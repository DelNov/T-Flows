!==============================================================================!
  program Generator
!------------------------------------------------------------------------------!
!   Block structured mesh generation and unstructured cell refinement.         !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Domain_Mod,  only: Domain_Type
  use Grid_Mod,    only: Grid_Type,                         &
                         Grid_Mod_Sort_Cells_Smart,         &
                         Grid_Mod_Sort_Faces_Smart,         &
                         Grid_Mod_Calculate_Wall_Distance,  &
                         Grid_Mod_Coarsen,                  &
                         Grid_Mod_Save_Cns,                 &
                         Grid_Mod_Save_Geo
  use Smooths_Mod, only: Smooths_Type
  use Refines_Mod, only: Refines_Type
!------------------------------------------------------------------------------!
  implicit none
!-----------------------------------[Locals]-----------------------------------!
  type(Domain_Type)  :: dom       ! domain to be used
  type(Grid_Type)    :: grid      ! grid which will be generated
  type(Smooths_Type) :: smooths   ! smoothing regions
  type(Refines_Type) :: refines   ! refinement regions and levels
  integer            :: c, s, n, lev
!==============================================================================!

  ! Open with a logo
  call Logo_Gen

  call Load_Domain               (dom, smooths, refines, grid)
  call Calculate_Node_Coordinates(dom, grid)
  call Distribute_Regions        (dom, grid)
  call Connect_Blocks            (dom, grid)
  call Connect_Periodicity       (dom, grid)
  call Connect_Copy              (dom, grid)

  ! From this point on, domain should not be used anymore
  call Determine_Grid_Connectivity(refines, grid, .false.)  ! trial run 
  call Calculate_Grid_Geometry    (grid, .false.)
  call Smooth_Grid                (smooths, grid)
  call Refine_Grid                (refines, grid)
  call Determine_Grid_Connectivity(refines, grid, .true.)   ! real run
  call Calculate_Grid_Geometry    (grid, .true.)

  call Grid_Mod_Sort_Cells_Smart       (grid)
  call Grid_Mod_Sort_Faces_Smart       (grid)
  call Grid_Mod_Calculate_Wall_Distance(grid)

  ! Prepare for saving
  do n = 1, grid % n_nodes
    grid % new_n(n) = n
  end do
  do c = -grid % n_bnd_cells, grid % n_cells
    grid % new_c(c) = c
  end do
  do s = 1, grid % n_faces
    grid % new_f(s) = s
  end do

  ! Coarsen the grid with METIS
  ! call Grid_Mod_Coarsen(grid, 1200, 4)

  !------------------------------!
  !   Save data for processing   !
  !------------------------------!
  call Grid_Mod_Save_Cns(grid, 0,             &
                         grid % n_nodes,      &
                         grid % n_cells,      &
                         grid % n_faces,      &
                         grid % n_bnd_cells,  &
                         0)

  call Grid_Mod_Save_Geo(grid, 0,             &
                         grid % n_faces,      &
                         0)

  !-----------------------------------------------------!
  !   Save grid for visualisation and post-processing   !
  !-----------------------------------------------------!

  ! Create output in vtu format
  call Save_Vtu_Cells(grid, 0,         &
                      grid % n_nodes,  &
                      grid % n_cells)
  call Save_Vtu_Faces(grid)

  ! Save all grid levels for visual inspection
  do lev = 1, grid % n_levels
    call Save_Vtu_Grid_Levels(grid, lev)
  end do

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
