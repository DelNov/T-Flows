!==============================================================================!
  subroutine Buoyancy_Forces(Flow, i)
!------------------------------------------------------------------------------!
!   Calculates buoyancy with or without Boussinesq hypothesis.  That means     !
!   it is used for both single phase flows with density assumed constant and   !
!   for multiphase flows with VOF and large density variations.                !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Field_Type), target :: Flow
  integer, intent(in)       :: i     ! component
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),   pointer :: Grid
  type(Var_Type),    pointer :: t
  integer                    :: c, c1, c2, s
  real                       :: xic1, xic2, temp_f, temp_ref_f
  real,              pointer :: grav_i
  real, contiguous,  pointer :: xic(:), xif(:), si(:), dxi(:)
  real, contiguous,  pointer :: cell_fi(:), face_fi(:)
  real, contiguous,  pointer :: dens_f(:)
!==============================================================================!

  call Work % Connect_Real_Face(dens_f)

  ! Take aliases
  Grid => Flow % pnt_grid
  t    => Flow % t

  ! Initialize
  ! Units: N/m^3
  if(i .eq. 1) then
    cell_fi => Flow % cell_fx;
    face_fi => Flow % face_fx;
    dxi     => Grid % dx;
    xic     => Grid % xc;
    xif     => Grid % xf;
    si      => Grid % sx;
    grav_i  => Flow % grav_x
  else if(i .eq. 2) then
    cell_fi => Flow % cell_fy;
    face_fi => Flow % face_fy;
    dxi     => Grid % dy;
    xic     => Grid % yc;
    xif     => Grid % yf;
    si      => Grid % sy;
    grav_i  => Flow % grav_y
  else if(i .eq. 3) then
    cell_fi => Flow % cell_fz;
    face_fi => Flow % face_fz;
    dxi     => Grid % dz;
    xic     => Grid % zc;
    xif     => Grid % zf;
    si      => Grid % sz;
    grav_i  => Flow % grav_z
  end if

  face_fi(:) = 0.0
  cell_fi(:) = 0.0

  !-------------------------------!
  !   For Boussinesq hypothesis   !
  !-------------------------------!
  if(Flow % buoyancy .eq. THERMALLY_DRIVEN) then

    ! If a per-cell reference temperature profile was requested
    ! (control file: REFERENCE_TEMPERATURE_PROFILE yes), capture it
    ! from the field once, on the first call - this is meant for cases
    ! with a non-trivial background stratification (e.g. penetrative
    ! convection), where buoyancy should be computed relative to the
    ! local ambient state rather than a single global constant.  This
    ! generalizes what used to be a hardcoded, Deardorff-specific
    ! formula here (t_ref = 21.5 + 45*z, matching that one case's
    ! initial linear profile) into something any stratified case can
    ! opt into just by setting the control-file flag.
    if(Flow % t_ref_profile_on .and. .not. allocated(Flow % t_ref_prof)) then
      allocate(Flow % t_ref_prof(-Grid % n_bnd_cells:Grid % n_cells))
      Flow % t_ref_prof(:) = t % n(:)
    end if

    do s = 1, Grid % n_faces
      c1  = Grid % faces_c(1, s)
      c2  = Grid % faces_c(2, s)

      ! Find temperature and density with linear interpolation.
      ! (At one point Hamo noticed that using straight "f" for
      !  interpolation gives too high buoyancy forces, and they
      !  have consequently been replaced by "fw" int. factors)
      temp_f =        Grid % fw(s)  * t % n(c1)  &
             + (1.0 - Grid % fw(s)) * t % n(c2)

      ! Density is constant for Boussinesq hypothesis
      ! this linear interpolation should do just fine
      dens_f(s) =        Grid % fw(s)  * Flow % density(c1)  &
                + (1.0 - Grid % fw(s)) * Flow % density(c2)

      ! Transform density into one assumed by Boussinesq
      ! Units: kg/m^3 * K * 1/K = kg/m^3
      if(Flow % t_ref_profile_on) then
        temp_ref_f =        Grid % fw(s)  * Flow % t_ref_prof(c1)  &
                   + (1.0 - Grid % fw(s)) * Flow % t_ref_prof(c2)
      else
        temp_ref_f = Flow % t_ref
      end if

      dens_f(s) = dens_f(s) * (temp_ref_f - temp_f) * Flow % beta

      ! Optional "sponge": zero buoyancy above a cutoff height, e.g.
      ! for a computational domain extended above the physically
      ! modelled region (control file: BUOYANCY_CUTOFF_HEIGHT).
      ! Defaults to HUGE, so this never triggers unless requested.
      if((Grid % zc(c1) + Grid % zc(c2)) * 0.5 > Flow % buoy_cutoff_z) then
        dens_f(s) = 0.0
      end if

    end do

  !---------------------------------------------------------------------!
  !   no Boussinesq hypothesis, this will include multiphase flow and   !
  !    probably very abbrupt changes in density, go for harmonic mean   !
  !---------------------------------------------------------------------!
  else

    ! Interpolate to all faces ...
    call Flow % Interpolate_To_Faces_Harmonic(dens_f, Flow % density)

    ! ... then substract the reference density from them all
    ! Note: dens_f(:) = dens_f(:) - Flow % dens_ref doesn't work for Intel
    do s = 1, Grid % n_faces
      dens_f(s) = dens_f(s) - Flow % dens_ref
    end do
  end if

  !--------------------!
  !   Buoyancy force   !
  !--------------------!

  ! Calculate buoyancy force over faces first (this is a no-brainer)
  ! Units: kg/m^3 * m/s^2 = kg/(m^2 s^2) = N/m^3
  do s = 1, Grid % n_faces
    face_fi(s) = dens_f(s) * grav_i
  end do

  ! Calculate
  do s = 1, Grid % n_faces
    c1   = Grid % faces_c(1, s)
    c2   = Grid % faces_c(2, s)
    xic1 = xic(c1)

    ! Units here: N/m^3 * m^3 = N
    ! (Don't worry yet, you will correct them in the end)
    cell_fi(c1) = cell_fi(c1) + face_fi(s) * si(s) * (xif(s) - xic1)

    if(c2 > 0) then
      xic2 = xic(c1) + dxi(s)

      ! Units here: N/m^3 * m^3 = N
      ! (Don't worry yet, you will correct them in the end)
      cell_fi(c2) = cell_fi(c2) - face_fi(s) * si(s) * (xif(s) - xic2)
    end if

  end do

  !---------------------------------------!
  !   Correct the units for body forces   !
  !---------------------------------------!
  do c = Cells_In_Domain()
    cell_fi(c) = cell_fi(c) / Grid % vol(c)
  end do

  call Grid % Exchange_Cells_Real(cell_fi)

  call Work % Disconnect_Real_Face(dens_f)

  end subroutine
