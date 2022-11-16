!==============================================================================!
  subroutine User_Mod_Beginning_Of_Simulation(Flow, Turb, Vof, Swarm,  &
                                              curr_dt, time)
!------------------------------------------------------------------------------!
!   This function is called at the beginning of simulation.                    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),    target :: Flow
  type(Turb_Type),     target :: Turb
  type(Vof_Type),      target :: Vof
  type(Swarm_Type),    target :: Swarm
  integer, intent(in)         :: curr_dt  ! time step
  real,    intent(in)         :: time     ! physical time
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: Grid
  real, contiguous, pointer :: por_real(:)   ! just for saving
  type(Forrest_Type)        :: forrest
  character(len=SL)         :: name_in
  integer                   :: fu, n_facets, v, c, f
  real                      :: cx, cy, cz, dot_prod
!==============================================================================!

  call Work % Connect_Real_Cell(por_real)

  ! Get alias
  Grid => Flow % pnt_grid

  ! Set the file name by hard coding
  name_in = 'forrest.stl'

  call File % Open_For_Reading_Ascii(name_in, fu)

  !--------------------------!
  !   Count all the facets   !
  !--------------------------!
  rewind(fu)
  n_facets = 0
  do
    call File % Read_Line(fu)
    if(Line % tokens(1) .eq. 'endsolid') exit
    if(Line % tokens(1) .eq. 'facet') n_facets = n_facets + 1
  end do
  print '(a38,i9)', '# Number of facets on the forrest:    ', n_facets

  allocate(forrest % facet(n_facets))

  !------------------------------------------------!
  !   Read all facet normals, vertex coordinates   !
  !   and compute facet centroids along the way    !
  !------------------------------------------------!
  rewind(fu)
  n_facets = 0
  do
    call File % Read_Line(fu)
    if(Line % tokens(1) .eq. 'endsolid') exit
    if(Line % tokens(1) .eq. 'facet') then
      n_facets = n_facets + 1
      read(Line % tokens(3), *) forrest % facet(n_facets) % nx
      read(Line % tokens(4), *) forrest % facet(n_facets) % ny
      read(Line % tokens(5), *) forrest % facet(n_facets) % nz
      call File % Read_Line(fu)                ! 'outer loop'
      forrest % facet(n_facets) % xc = 0.0
      forrest % facet(n_facets) % yc = 0.0
      forrest % facet(n_facets) % zc = 0.0
      do v = 1, 3
        call File % Read_Line(fu)              ! 'vertex 1, 2 and 3'
        read(Line % tokens(2), *) forrest % facet(n_facets) % x(v)
        read(Line % tokens(3), *) forrest % facet(n_facets) % y(v)
        read(Line % tokens(4), *) forrest % facet(n_facets) % z(v)
        forrest % facet(n_facets) % xc = forrest % facet(n_facets) % xc    &
                                  + forrest % facet(n_facets) % x(v)
        forrest % facet(n_facets) % yc = forrest % facet(n_facets) % yc    &
                                  + forrest % facet(n_facets) % y(v)
        forrest % facet(n_facets) % zc = forrest % facet(n_facets) % zc    &
                                  + forrest % facet(n_facets) % z(v)
      end do
      forrest % facet(n_facets) % xc = forrest % facet(n_facets) % xc * ONE_THIRD
      forrest % facet(n_facets) % yc = forrest % facet(n_facets) % yc * ONE_THIRD
      forrest % facet(n_facets) % zc = forrest % facet(n_facets) % zc * ONE_THIRD
      call File % Read_Line(fu)                ! 'endloop'
    end if
  end do
  print *, '# Read all forrest facets!'

  !-----------------------------------------!
  !   Find which cells are in the forrest   !
  !-----------------------------------------!
  allocate(porosity(Grid % n_cells));  porosity(:) = .false.

  do c = 1, Grid % n_cells

    ! Assume cell is in the forrest
    porosity(c) = .true.
    do f = 1, n_facets

      ! Vector connecting facet centroid with cell centroid
      cx = forrest % facet(f) % xc - Grid % xc(c)
      cy = forrest % facet(f) % yc - Grid % yc(c)
      cz = forrest % facet(f) % zc - Grid % zc(c)

      dot_prod = cx * forrest % facet(f) % nx  &
               + cy * forrest % facet(f) % ny  &
               + cz * forrest % facet(f) % nz

      ! First time this is negative, cell is not in the forrest
      if(dot_prod < 0) then
        porosity(c) = .false.
        goto 1
      end if
    end do
1   continue

  end do

  por_real(1:Grid % n_cells) = 0.0
  do c = 1, Grid % n_cells
    if(porosity(c)) por_real(c) = 1.0;
  end do

  call Grid % Save_Debug_Vtu('porosity',             &
                             scalar_cell=por_real,   &
                             scalar_name='porosity')

  call Work % Disconnect_Real_Cell(por_real)

  end subroutine
