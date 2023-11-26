!==============================================================================!
  pure integer function N_Procs()
!------------------------------------------------------------------------------!
!>  A straightforward utility that returns the total number of processors
!>  involved in the current execution. It simply retrieves and returns the
!>  value of private n_processors from the global Comm_Type object Global.
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  N_Procs = Global % n_processors

  end function
