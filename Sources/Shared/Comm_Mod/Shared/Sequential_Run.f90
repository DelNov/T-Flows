!==============================================================================!
  pure logical function Sequential_Run()
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  Sequential_Run = Communicator % n_processors < 2

  end function