!==============================================================================!
  pure integer function N_Procs()
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  N_Procs = Communicator % n_processors

  end function