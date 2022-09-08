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
  type(Grid_Type), pointer :: Grid
  integer                  :: c1, c2, s, i_nod, j_nod
  real                     :: xf(MAX_FACES_N_NODES),  &
                              yf(MAX_FACES_N_NODES),  &
                              zf(MAX_FACES_N_NODES)
  real                     :: x_cell_1, y_cell_1, z_cell_1
  real                     :: x_cell_2, y_cell_2, z_cell_2
  real                     :: ix, iy, iz, ixy, ixz, iyz
!==============================================================================!

  ! Make sure you compiled it
  print *, '# In User_Mod_Beginning_Of_Simulation'

  ! Take alias
  Grid => Flow % pnt_grid

  ! Allocate memory for user arrays
  allocate(cell_inertia % x (-Grid % n_bnd_cells : Grid % n_cells))
  allocate(cell_inertia % y (-Grid % n_bnd_cells : Grid % n_cells))
  allocate(cell_inertia % z (-Grid % n_bnd_cells : Grid % n_cells))
  allocate(cell_inertia % xy(-Grid % n_bnd_cells : Grid % n_cells))
  allocate(cell_inertia % xz(-Grid % n_bnd_cells : Grid % n_cells))
  allocate(cell_inertia % yz(-Grid % n_bnd_cells : Grid % n_cells))
  print *, '# Allocated memory for cell inertia tensor components'

  ! Initialize cell inertias
  cell_inertia % x (:) = 0.0
  cell_inertia % y (:) = 0.0
  cell_inertia % z (:) = 0.0
  cell_inertia % xy(:) = 0.0
  cell_inertia % xz(:) = 0.0
  cell_inertia % yz(:) = 0.0

  do s = 1, Grid % n_faces
    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)

    do i_nod = 1, Grid % faces_n_nodes(s)  ! for all face types
      xf(i_nod) = Grid % xn(Grid % faces_n(i_nod,s))
      yf(i_nod) = Grid % yn(Grid % faces_n(i_nod,s))
      zf(i_nod) = Grid % zn(Grid % faces_n(i_nod,s))
    end do

    ! First cell
    x_cell_1 = Grid % xc(c1)
    y_cell_1 = Grid % yc(c1)
    z_cell_1 = Grid % zc(c1)

    ! Browse through faces's nodes and add tetrahedron by tetrahedron
    do i_nod = 1, Grid % faces_n_nodes(s)
      j_nod = i_nod + 1;  if(j_nod > Grid % faces_n_nodes(s)) j_nod = 1

      call Math % Tet_Inertia(Grid % xf(s), Grid % yf(s), Grid % zf(s),  &
                              xf(i_nod),    yf(i_nod),    zf(i_nod),     &
                              xf(j_nod),    yf(j_nod),    zf(j_nod),     &
                              x_cell_1,     y_cell_1,     z_cell_1,      &
                              ix, iy, iz, ixy, ixz, iyz,                 &
                              around_node = 4)

      ! Update inertia in the first cell (not complete)
      cell_inertia % x (c1) = cell_inertia % x (c1) + ix
      cell_inertia % y (c1) = cell_inertia % y (c1) + iy
      cell_inertia % z (c1) = cell_inertia % z (c1) + iz
      cell_inertia % xy(c1) = cell_inertia % xy(c1) + ixy
      cell_inertia % xz(c1) = cell_inertia % xz(c1) + ixz
      cell_inertia % yz(c1) = cell_inertia % yz(c1) + iyz

    end do  ! i_nod

    ! Second cell
    if(c2 > 0) then
      x_cell_2 = Grid % xc(c1) + Grid % dx(s)
      y_cell_2 = Grid % yc(c1) + Grid % dy(s)
      z_cell_2 = Grid % zc(c1) + Grid % dz(s)

      ! Browse through faces's nodes and add tetrahedron by tetrahedron
      do i_nod = 1, Grid % faces_n_nodes(s)
        j_nod = i_nod + 1;  if(j_nod > Grid % faces_n_nodes(s)) j_nod = 1

        call Math % Tet_Inertia(Grid % xf(s), Grid % yf(s), Grid % zf(s),  &
                                xf(j_nod),    yf(j_nod),    zf(j_nod),     &
                                xf(i_nod),    yf(i_nod),    zf(i_nod),     &
                                x_cell_2,     y_cell_2,     z_cell_2,      &
                                ix, iy, iz, ixy, ixz, iyz,                 &
                                around_node = 4)

        ! Update inertia in the second cell (not complete, continue)
        cell_inertia % x (c2) = cell_inertia % x (c2) + ix
        cell_inertia % y (c2) = cell_inertia % y (c2) + iy
        cell_inertia % z (c2) = cell_inertia % z (c2) + iz
        cell_inertia % xy(c2) = cell_inertia % xy(c2) + ixy
        cell_inertia % xz(c2) = cell_inertia % xz(c2) + ixz
        cell_inertia % yz(c2) = cell_inertia % yz(c2) + iyz

      end do  ! i_nod
    end if

  end do
  print *, '# Cell inertia calculated !'

  ! Refresh buffers for M % sav before discretizing for pressure
  call Grid % Exchange_Cells_Real(cell_inertia % x )
  call Grid % Exchange_Cells_Real(cell_inertia % y )
  call Grid % Exchange_Cells_Real(cell_inertia % z )
  call Grid % Exchange_Cells_Real(cell_inertia % xy)
  call Grid % Exchange_Cells_Real(cell_inertia % xz)
  call Grid % Exchange_Cells_Real(cell_inertia % yz)

  !--------------------------------------------------!
  !   Check the computed cells' moments of inertia   !
  !--------------------------------------------------!
  call Grid % Save_Debug_Vtu('cell-inertia',                       &
                             tensor_cell = (/cell_inertia % x,     &
                                             cell_inertia % y,     &
                                             cell_inertia % z,     &
                                             cell_inertia % xy,    &
                                             cell_inertia % xz,    &
                                             cell_inertia % yz/),  &
                             tensor_name = 'Cell Inertia')

  !----------------------------------------------------------------!
  !   Check the computed components of cells' moments of inertia   !
  !----------------------------------------------------------------!
  ! call Grid % Save_Debug_Vtu('ix',                            &
  !                            scalar_cell = cell_inertia % x,  &
  !                            scalar_name = 'Cell Inertia X')
  ! call Grid % Save_Debug_Vtu('iy',                            &
  !                            scalar_cell = cell_inertia % y,  &
  !                            scalar_name = 'Cell Inertia Y')
  ! call Grid % Save_Debug_Vtu('iz',                            &
  !                            scalar_cell = cell_inertia % z,  &
  !                            scalar_name = 'Cell Inertia Z')

  call Comm_Mod_End
  stop

  end subroutine
