#include "../Shared/Assert.h90"

!==============================================================================!
  program Generate_Prog
!------------------------------------------------------------------------------!
!>  Generate's main function.  Generate is a part of T-Flows and is designed
!>  for block structured mesh generation and unstructured cell refinement.
!------------------------------------------------------------------------------!
!   Program Flow:                                                              !
!                                                                              !
!   * Initialization: The program begins with a profiling start, followed by   !
!     displaying its logo.                                                     !
!   * Input reading: It reads the input file using Generate % Load_Dom, which  !
!     loads domain specifications, refinement criteria, and smoothing          !
!     parameters.                                                              !
!   * Domain handling: The domain's node coordinates are calculated, and       !
!     ranges are distributed. Blocks are connected, and periodicity in the     !
!     domain is established.                                                   !
!   * Refinement and smoothing: The grid undergoes a trial connectivity run,   !
!     followed by geometry calculation, smoothing, cell marking for refinement,!
!      and a final connectivity pass.                                          !
!   * Grid processing: The program sorts cells and faces, calculates wall      !
!     distances, and associates cells with faces.                              !
!   * Data saving and visualization: The processed grid data is saved in .cfn  !
!     and .dim formats for further processing with T-Flows and in VTU format   !
!     for visualization. A template control file is also created.              !
!   * Additional analysis: The program checks for grid homogeneity and saves   !
!     1D probes, which are helpful for channel flow analyses.                  !
!   * Final output: Generate statistics are printed, and profiling data are    !
!     finalized.                                                               !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Generate_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Interfaces]---------------------------------!
  interface
    include '../Shared/Probe_1d_Nodes.h90'
  end interface
!-----------------------------------[Locals]-----------------------------------!
  type(Domain_Type)  :: Dom       ! domain to be used
  type(Grid_Type)    :: Grid      ! Grid which will be generated
  type(Smooths_Type) :: smooths   ! smoothing regions
  type(Refines_Type) :: refines   ! refinement regions and levels
  integer            :: c, s      ! cell counter
!==============================================================================!

  ! Initialize program profler
  call Profiler % Start('Main')

  ! Open with a logo
  call Generate % Logo_Gen()

  !-------------------------!
  !   Read the input file   !
  !-------------------------!
  call Generate % Load_Dom(Dom, smooths, refines, Grid)

  !-----------------------!
  !   Handle the domain   !
  !-----------------------!
  call Dom % Calculate_Node_Coordinates(Grid)
  call Dom % Distribute_Ranges         (Grid)
  call Dom % Connect_Blocks            (Grid)
  call Dom % Connect_Periodicity       (Grid)

  !--------------------------------!
  !   From this point on, domain   !
  !   should not be used anymore   !
  !--------------------------------!
  call Refines_Mod_Connectivity     (refines, Grid, .false.)  ! trial run
  call Generate % Calculate_Geometry(Grid,          .false.)
  call Smooths_Mod_Grid             (smooths, Grid)
  call Refines_Mod_Mark_Cells       (refines, Grid)
  call Refines_Mod_Connectivity     (refines, Grid, .true.)   ! real run
  call Generate % Calculate_Geometry(Grid,          .true.)

  call Grid % Sort_Cells_By_Coordinates()
  call Grid % Sort_Faces_By_Region     ()
  call Grid % Calculate_Wall_Distance  ()
  call Grid % Find_Cells_Faces         ()

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
  do s = 1, Grid % n_faces + Grid % n_shadows
    if(Grid % faces_s(s) .ne. 0) then
      Assert(Grid % faces_s(Grid % faces_s(s)) .eq. s)
    end if
  end do
  ! Similar note is in Convert, also called Note #1

  !------------------------------!
  !   Save data for processing   !
  !------------------------------!
  call Grid % Save_Cfn((/0, 0/),            &
                       Grid % n_nodes,      &
                       Grid % n_cells,      &
                       Grid % n_faces,      &
                       Grid % n_shadows,    &
                       Grid % n_bnd_cells)
  call Grid % Save_Dim((/0, 0/))

  !-----------------------------------------------------!
  !   Save Grid for visualisation and post-processing   !
  !-----------------------------------------------------!

  ! Create output in vtu format
  call Grid % Save_Vtu_Cells((/0, 0/),        &
                             Grid % n_nodes,  &
                             Grid % n_cells)
  call Grid % Save_Vtu_Faces((/0, 0/))
  call Grid % Save_Vtu_Faces((/0, 0/), plot_inside=.false.)
  call Grid % Save_Vtu_Faces((/0, 0/), plot_shadows=.true.)

  ! Create a template control file for this domain
  call Grid % Write_Template_Control_File()

  ! Check homogeneity of the grid in three coordinate directions
  ! Here, in generator, there is no sense to enforce uniformity
  call Grid % Search_Coordinate_Clusters(nodal           = .true.,   &
                                         enforce_uniform = .false.)

  ! Save the 1d probe (good for the channel flow)
  call Probe_1d_Nodes(Grid)

  ! Write something on the screen
  call Generate % Print_Generate_Statistics(Grid)

  ! Finalize program profler
  call Profiler % Stop('Main')
  call Profiler % Statistics(indent=1)

  end program
