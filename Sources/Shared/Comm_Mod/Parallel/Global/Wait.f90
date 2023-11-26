!==============================================================================!
  subroutine Wait(Global)
!------------------------------------------------------------------------------!
!>  This subroutine acts as a synchronization point in parallel computing.
!>  It ensures that all processes in a parallel computation reach this point
!>  before any of them proceedes further, thus effectivelly preventing the
!>  "race conditions".
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type) :: Global  !! global communication object
!-----------------------------------[Locals]-----------------------------------!
  integer :: error
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Global)
!==============================================================================!

  call Mpi_Barrier(MPI_COMM_WORLD, error)

  end subroutine
