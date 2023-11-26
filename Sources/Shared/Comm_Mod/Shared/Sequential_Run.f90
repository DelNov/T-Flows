!==============================================================================!
  pure logical function Sequential_Run()
!------------------------------------------------------------------------------!
!>  The Sequential_Run function is designed to determine the execution mode of
!>  a program. It returns false if program is run in parallel (on multiple
!>  processors) and true if it runs sequentially (on one processor).
!------------------------------------------------------------------------------!
!   Note that This function returns:                                           !
!                                                                              !
!   * .true.  - for sequential programs: Generate, Convert and Divide          !
!   * .true.  - for parallel program Process, if compiled without MPI          !
!   * .true.  - for parallel program Process, if ran on one processor          !
!   * .false. - for parallel program Process, if ran on more processors        !
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  Sequential_Run = Global % n_processors < 2

  end function
