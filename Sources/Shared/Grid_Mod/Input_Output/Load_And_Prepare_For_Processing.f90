!==============================================================================!
  subroutine Load_And_Prepare_For_Processing(Grid, d)
!------------------------------------------------------------------------------!
!>  This subroutine is designed for use in the Process program and performs
!>  several critical steps to prepare the computational grid for processing. It
!>  includes loading grid data, setting up OpenMP threads, and establishing
!>  communication patterns essential for parallel computations
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * Checking program context:                                                !
!     - Asserts that the subroutine is called from the Process program to      !
!       ensure it's used in the correct context.                               !
!   * Load grid data:                                                          !
!     - Reads the problem name corresponding to the domain 'd'.                !
!     - Loads the finite volume grid data using the Load_Cfn and Load_Dim      !
!       subroutines, which handle connectivity and geometrical aspects of the  !
!       grid respectively.                                                     !
!   * OpenMP thread setup:                                                     !
!     - Reads the maximum number of threads from the control settings and      !
!       configures OpenMP threads accordingly using Determine_Threads.         !
!   * Calculate face geometry:                                                 !
!     - Computes geometrical properties of the grid faces, which are crucial   !
!       for various numerical calculations in the simulation.                  !
!   * Establish communication patterns:                                        !
!     - Sets up communication patterns for cells, which is vital for parallel  !
!       processing and data exchange among different processors.               !
!     - Prepares maps for backup purposes, ensuring robustness in data         !
!       handling across different computational nodes.                         !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type), intent(inout) :: Grid  !! computational grid
  integer,          intent(in)    :: d     !! domain rank (number)
!==============================================================================!

  Assert(PROGRAM_NAME(1:7) .eq. 'Process')

  call Control % Read_Problem_Name(problem_name(d))

  ! Load the finite volume Grid
  call Grid % Load_Cfn((/This_Proc(), N_Procs()/), domain=d)
  call Grid % Load_Dim((/This_Proc(), N_Procs()/), domain=d)

  ! Determine threads for OpenMP runs
  call Control % Max_Threads(Grid % Omp % d_threads, .true.)
  call Grid % Determine_Threads()

  call Grid % Calculate_Face_Geometry()

  ! Find communication patterns
  call Grid % Form_Cells_Comm()
  call Grid % Form_Maps_For_Backup()

  end subroutine

