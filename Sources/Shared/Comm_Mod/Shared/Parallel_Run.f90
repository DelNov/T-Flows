!==============================================================================!
  pure logical function Parallel_Run()
!------------------------------------------------------------------------------!
!>  The Parallel_Run function is designed to determine the execution mode of
!>  a program. It returns true if program is run in parallel (on multiple
!>  processors) and false if it runs sequentially (on one processor).
!------------------------------------------------------------------------------!
!   Note that This function returns:                                           !
!                                                                              !
!   * .false. - for sequential programs: Generate, Convert and Divide          !
!   * .false. - for parallel program Process, if compiled without MPI          !
!   * .false. - for parallel program Process, if ran on one processor          !
!   * .true.  - for parallel program Process, if ran on more processors        !
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  Parallel_Run = Global % n_processors > 1

  end function
