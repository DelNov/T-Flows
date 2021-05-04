!==============================================================================!
  subroutine Field_Mod_Buoyancy_Forces(flow, i)
!------------------------------------------------------------------------------!
!   Calculates buoyancy with or without Boussinesq hypothesis.  That means     !
!   it is used for both single phase flows with density assumed constant and   !
!   for multiphase flows with VOF and large density variations.                !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Work_Mod, only: dens_f => r_face_01
!------------------------------------------------------------------------------!
!   When using Work_Mod, calling sequence should be outlined                   !
!                                                                              !
!   Main_Pro                                    (allocates Work_Mod)           !
!     |                                                                        !
!     +----> Compute_Momentum                   (doesn't use Work_Mod)         !
!             |                                                                !
!             +----> Field_Mod_Buoyancy_Forces  (safe to use r_face_01)        !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type), target :: flow
  integer, intent(in)      :: i     ! component
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),   pointer :: grid
  type(Var_Type),    pointer :: t
  integer                    :: c, c1, c2, s
  real                       :: xic1, xic2, temp_f
  real,              pointer :: grav_i
  real, contiguous,  pointer :: xic(:), xif(:), si(:), dxi(:)
  real, contiguous,  pointer :: cell_fi(:), face_fi(:)
!==============================================================================!

  ! Take aliases
  grid => flow % pnt_grid
  t    => flow % t

  ! Initialize
  ! Units: N/m^3
  if(i .eq. 1) then
    cell_fi => flow % cell_fx;
    face_fi => flow % face_fx;
    dxi     => grid % dx;
    xic     => grid % xc;
    xif     => grid % xf;
    si      => grid % sx;
    grav_i  => grav_x
  else if(i .eq. 2) then
    cell_fi => flow % cell_fy;
    face_fi => flow % face_fy;
    dxi     => grid % dy;
    xic     => grid % yc;
    xif     => grid % yf;
    si      => grid % sy;
    grav_i  => grav_y
  else if(i .eq. 3) then
    cell_fi => flow % cell_fz;
    face_fi => flow % face_fz;
    dxi     => grid % dz;
    xic     => grid % zc;
    xif     => grid % zf;
    si      => grid % sz;
    grav_i  => grav_z
  end if

  face_fi(:) = 0.0
  cell_fi(:) = 0.0

  !-------------------------------!
  !   For Boussinesq hypothesis   !
  !-------------------------------!
  if(boussinesq) then
    do s = 1, grid % n_faces
      c1  = grid % faces_c(1, s)
      c2  = grid % faces_c(2, s)

      ! Find temperature and density with linear interpolation.
      ! (One could consider a more accurate interpolation with
      !  gradients, or even just running it through Grad_Gauss
      !  to get iterativelly improved face temperature values)
      temp_f =        grid % f(s)  * t % n(c1)  &
             + (1.0 - grid % f(s)) * t % n(c2)

      ! Density is constant for Boussinesq hypothesis
      ! this linear interpolation should do just fine
      dens_f(s) =        grid % fw(s)  * flow % density(c1)  &
                + (1.0 - grid % fw(s)) * flow % density(c2)

      ! Transform density into one assumed by Boussinesq
      ! Units: kg/m^3 * K * 1/K = kg/m^3
      dens_f(s) = dens_f(s) * (flow % t_ref - temp_f) * flow % beta
    end do

  !---------------------------------------------------------------------!
  !   no Boussinesq hypothesis, this will include multiphase flow and   !
  !    probably very abbrupt changes in density, go for harmonic mean   !
  !---------------------------------------------------------------------!
  else

    ! Interpolate to all faces ...
    call Field_Mod_Interpolate_To_Faces_Harmonic(flow, dens_f, flow % density)

    ! ... then substract the reference density from them all
    dens_f(:) = dens_f(:) - flow % dens_ref
  end if

  !--------------------!
  !   Buoyancy force   !
  !--------------------!

  ! Calculate
  do s = 1, grid % n_faces
    c1   = grid % faces_c(1, s)
    c2   = grid % faces_c(2, s)
    xic1 = xic(c1)

    ! Units in next three expressions which follow
    ! kg/m^3 * m/s^2 * m^2 * m = kg m/s^2 = N
    ! (Don't worry yet, you will correct them in the end)
    cell_fi(c1) = cell_fi(c1)                                  &
                + dens_f(s) * grav_i * si(s) * (xif(s) - xic1)

    if(c2 > 0) then
      xic2 = xic(c1) + dxi(s)

      ! Units in next three expressions which follow
      ! kg/m^3 * m/s^2 * m^2 * m = kg m/s^2 = N
      ! (Don't worry yet, you will correct them in the end)
      cell_fi(c2) = cell_fi(c2)                                &
                  - dens_f(s) * grav_i * si(s) * (xif(s) - xic2)
    end if

    ! Units: kg/m^3 * m/s^2 = kg/(m^2 s^2) = N/m^3
    face_fi(s) = face_fi(s) + dens_f(s) * grav_i

  end do

  !---------------------------------------!
  !   Correct the units for body forces   !
  !---------------------------------------!
  do c = 1, grid % n_cells
    cell_fi(c) = cell_fi(c) / grid % vol(c)
  end do

  call Grid_Mod_Exchange_Cells_Real(grid, cell_fi)

  end subroutine
