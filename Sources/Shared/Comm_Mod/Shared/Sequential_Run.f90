!==============================================================================!
  pure logical function Sequential_Run()
!------------------------------------------------------------------------------!
!   This function returns:                                                     !
!     .true.  - for sequential programs: Generate, Convert and Divide          !
!     .true.  - for parallel program Process, if compiled without MPI          !
!     .true.  - for parallel program Process, if ran on one processor          !
!     .false. - for parallel program Process, if ran on more processors        !
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  Sequential_Run = Global % n_processors < 2

  end function
