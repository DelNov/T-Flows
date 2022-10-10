!==============================================================================!
  program Generator
!------------------------------------------------------------------------------!
!   Block structured mesh generation and unstructured cell refinement.         !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Domain_Mod
  use Smooths_Mod
  use Refines_Mod
  use Swap_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Interfaces]---------------------------------!
  interface
    include 'Calculate_Geometry.h90'
    include 'Load_Dom.h90'
    include 'Logo_Gen.h90'
    include '../Shared/Probe_1d_Nodes.h90'
    include '../Shared/Probe_2d.h90'
  end interface
!-----------------------------------[Locals]-----------------------------------!
  type(Domain_Type)  :: dom       ! domain to be used
  type(Grid_Type)    :: Grid      ! Grid which will be generated
  type(Smooths_Type) :: smooths   ! smoothing regions
  type(Refines_Type) :: refines   ! refinement regions and levels
  integer            :: c         ! cell counter
!==============================================================================!

  ! Open with a logo
  call Logo_Gen

  !-------------------------!
  !   Read the input file   !
  !-------------------------!
  call Load_Dom(dom, smooths, refines, Grid)

  !-----------------------!
  !   Handle the domain   !
  !-----------------------!
  call Domain_Mod_Calculate_Node_Coordinates(dom, Grid)
  call Domain_Mod_Distribute_Regions        (dom, Grid)
  call Domain_Mod_Connect_Blocks            (dom, Grid)
  call Domain_Mod_Connect_Periodicity       (dom, Grid)

  !--------------------------------!
  !   From this point on, domain   !
  !   should not be used anymore   !
  !--------------------------------!
  call Refines_Mod_Connectivity(refines, Grid, .false.)  ! trial run
  call Calculate_Geometry      (Grid,          .false.)
  call Smooths_Mod_Grid        (smooths, Grid)
  call Refines_Mod_Mark_Cells  (refines, Grid)
  call Refines_Mod_Connectivity(refines, Grid, .true.)   ! real run
  call Calculate_Geometry      (Grid,          .true.)

  call Grid % Sort_Cells_Smart       ()
  call Grid % Sort_Faces_Smart       ()
  call Grid % Calculate_Wall_Distance()
  call Grid % Find_Cells_Faces       ()

  ! Prepare for saving
  call Grid % Initialize_New_Numbers()

  ! Make cell numberig compatible with VTU format
  do c = 1, Grid % n_cells
    call Swap_Int(Grid % cells_n(3,c), Grid % cells_n(4,c))
    call Swap_Int(Grid % cells_n(7,c), Grid % cells_n(8,c))
  end do

  ! Note #1 about shadows:
  ! At this point you have Grid % n_faces faces and Grid % n_shadows (on top)
  ! and they are pointing to each other.  Besides, both real face and its
  ! shadow have the same c1 and c2, both inside cells with positive indices
  ! Real faces which do not have shadows have "0" for shadow.
  ! Checked like this: do s = 1, Grid % n_faces + Grid % n_shadows
  ! Checked like this:   write(20, '(99i9)') s, Grid % faces_s(s)
  ! Checked like this: end do
  ! Similar note is in Convert, also called Note #1

  !------------------------------!
  !   Save data for processing   !
  !------------------------------!
  call Grid % Save_Cfn(0,                   &
                       Grid % n_nodes,      &
                       Grid % n_cells,      &
                       Grid % n_faces,      &
                       Grid % n_shadows,    &
                       Grid % n_bnd_cells)
  call Grid % Save_Dim(0)

  !-----------------------------------------------------!
  !   Save Grid for visualisation and post-processing   !
  !-----------------------------------------------------!

  ! Create output in vtu format
  call Grid % Save_Vtu_Cells(0,               &
                             Grid % n_nodes,  &
                             Grid % n_cells)
  call Grid % Save_Vtu_Faces()
  call Grid % Save_Vtu_Faces(plot_shadows=.true.)

  ! Create a template control file for this domain
  call Grid % Write_Template_Control_File()

  ! Save the 1d probe (good for the channel flow)
  call Probe_1d_Nodes(Grid)

  ! Save the 2d probe
  call Probe_2d(Grid)

  ! Write something on the screen
  call Print_Grid_Statistics(Grid)

  end program
