!==============================================================================!
  program Generator
!------------------------------------------------------------------------------!
!   Block structured mesh generation and unstructured cell refinement.         !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Gen_Mod
  use Domain_Mod  ! domain as defined in ".dom" file.
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!-----------------------------------[Locals]-----------------------------------!
  type(Domain_Type) :: dom    ! domain to be used
  type(Grid_Type)   :: grid   ! grid which will be generated
  integer           :: c, s, n
!==============================================================================!

  ! Open with a logo
  call Logo_Gen

  call Load_Domain               (dom, grid)
  call Calculate_Node_Coordinates(dom, grid)
  call Distribute_Regions        (dom, grid)
  call Connect_Blocks            (dom, grid)
  call Connect_Periodicity       (dom, grid)
  call Connect_Copy              (dom, grid)

  ! From this point on, domain is not used anymore
  call Determine_Grid_Connectivity(grid, .false.)  ! trial run 
  call Calculate_Grid_Geometry    (grid, .false.)
  call Smooth_Grid                (grid)
  call Refine_Grid                (grid)
  call Determine_Grid_Connectivity(grid, .true.)   ! real run
  call Calculate_Grid_Geometry    (grid, .true.)

  ! Prepare for saving
  do n = 1, grid % n_nodes
    new_n(n) = n
  end do
  do c = -grid % n_bnd_cells, grid % n_cells
    new_c(c) = c
  end do
  do s = 1, grid % n_faces
    new_f(s) = s
  end do

  ! Coarsen the grid with METIS
  ! Not implemented yet: call Grid_Mod_Coarsen(grid)

  !------------------------------!
  !   Save data for processing   !
  !------------------------------!
  call Save_Cns_Geo(grid, 0,             &
                    grid % n_nodes,      &
                    grid % n_cells,      &
                    grid % n_faces,      &
                    grid % n_bnd_cells,  &
                    0, 0)  ! saved data for processing

  !-----------------------------------------------------!
  !   Save grid for visualisation and post-processing   !
  !-----------------------------------------------------!

  ! Create output in vtu format
  call Save_Vtu_Cells(grid, 0,         &
                      grid % n_nodes,  &
                      grid % n_cells)
  call Save_Vtu_Faces(grid)

  ! Save links for checking
  call Save_Vtu_Links(grid, 0,             &
                      grid % n_nodes,      &
                      grid % n_cells,      &
                      grid % n_faces,      &
                      grid % n_bnd_cells,  &
                      0)

  ! Try to save in CGNS format, it might work
  call Save_Cgns_Cells(grid, 0) 

  ! Save the 1d probe (good for the channel flow)
  call Probe_1d_Nodes(grid)

  ! Save the 2d probe
  call Probe_2d(grid)

  ! Write something on the screen
  call Print_Grid_Statistics(grid)

  end program
