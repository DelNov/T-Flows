!==============================================================================!
  program Generator
!------------------------------------------------------------------------------!
!   Block structured mesh generation and unstructured cell refinement.         !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Domain_Mod
  use Smooths_Mod
  use Refines_Mod
  use Save_Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!-----------------------------------[Locals]-----------------------------------!
  type(Domain_Type)  :: dom       ! domain to be used
  type(Grid_Type)    :: grid      ! grid which will be generated
  type(Smooths_Type) :: smooths   ! smoothing regions
  type(Refines_Type) :: refines   ! refinement regions and levels
  integer            :: c         ! cell counter
!==============================================================================!

  ! Open with a logo
  call Logo_Gen

  !-------------------------!
  !   Read the input file   !
  !-------------------------!
  call Load_Dom(dom, smooths, refines, grid)

  !-----------------------!
  !   Handle the domain   !
  !-----------------------!
  call Domain_Mod_Calculate_Node_Coordinates(dom, grid)
  call Domain_Mod_Distribute_Regions        (dom, grid)
  call Domain_Mod_Connect_Blocks            (dom, grid)
  call Domain_Mod_Connect_Periodicity       (dom, grid)

  !--------------------------------!
  !   From this point on, domain   !
  !   should not be used anymore   !
  !--------------------------------!
  call Refines_Mod_Connectivity(refines, grid, .false.)  ! trial run
  call Calculate_Geometry      (grid,          .false.)
  call Smooths_Mod_Grid        (smooths, grid)
  call Refines_Mod_Mark_Cells  (refines, grid)
  call Refines_Mod_Connectivity(refines, grid, .true.)   ! real run
  call Calculate_Geometry      (grid,          .true.)

  call Grid_Mod_Sort_Cells_Smart       (grid)
  call Grid_Mod_Sort_Faces_Smart       (grid)
  call Grid_Mod_Calculate_Wall_Distance(grid)
  call Grid_Mod_Find_Cells_Faces       (grid)

  ! Prepare for saving
  call Grid_Mod_Initialize_New_Numbers(grid)

  ! Make cell numberig compatible with VTU format
  do c = 1, grid % n_cells
    call Swap_Int(grid % cells_n(3,c), grid % cells_n(4,c))
    call Swap_Int(grid % cells_n(7,c), grid % cells_n(8,c))
  end do

  ! Note #1 about shadows:
  ! At this point you have grid % n_faces faces and grid % n_shadows (on top)
  ! and they are pointing to each other.  Besides, both real face and its
  ! shadow have the same c1 and c2, both inside cells with positive indices
  ! Real faces which do not have shadows have "0" for shadow.
  ! Checked like this: do s = 1, grid % n_faces + grid % n_shadows
  ! Checked like this:   write(20, '(99i9)') s, grid % faces_s(s)
  ! Checked like this: end do
  ! Similar note is in Convert, also called Note #1

  !------------------------------!
  !   Save data for processing   !
  !------------------------------!
  call Grid_Mod_Save_Cfn(grid, 0,             &
                         grid % n_nodes,      &
                         grid % n_cells,      &
                         grid % n_faces,      &
                         grid % n_shadows,    &
                         grid % n_bnd_cells)

  call Grid_Mod_Save_Dim(grid, 0)

  !-----------------------------------------------------!
  !   Save grid for visualisation and post-processing   !
  !-----------------------------------------------------!

  ! Create output in vtu format
  call Save_Vtu_Cells(grid, 0,         &
                      grid % n_nodes,  &
                      grid % n_cells)
  call Save_Vtu_Faces(grid)
  call Save_Vtu_Faces(grid, plot_shadows=.true.)

  ! Try to save in CGNS format, it might work
  call Save_Cgns_Cells(grid, 0) 

  ! Save the 1d probe (good for the channel flow)
  call Probe_1d_Nodes(grid)

  ! Save the 2d probe
  call Probe_2d(grid)

  ! Write something on the screen
  call Print_Grid_Statistics(grid)

  end program
