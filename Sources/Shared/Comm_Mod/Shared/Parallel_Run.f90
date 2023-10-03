!==============================================================================!
  pure logical function Parallel_Run()
!------------------------------------------------------------------------------!
!   This function returns:                                                     !
!     .false. - for sequential programs: Generate, Convert and Divide          !
!     .false. - for parallel program Process, if compiled without MPI          !
!     .false. - for parallel program Process, if ran on one processor          !
!     .true.  - for parallel program Process, if ran on more processors        !
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  Parallel_Run = Global % n_processors > 1

  end function
