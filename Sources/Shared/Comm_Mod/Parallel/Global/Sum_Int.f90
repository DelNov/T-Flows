!==============================================================================!
  subroutine Sum_Int(Global, phi)
!------------------------------------------------------------------------------!
!>  The subroutine Sum_Int in a parallel computing framework calculates
!>  the global sum of integer value across all processors. This subroutine is
!>  vital in scenarios where an aggregate sum of a distributed integer value is
!>  required for collective decision-making or further computational processes.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type), intent(in)    :: Global  !! global communication object
  integer,          intent(inout) :: phi     !! global sum
!-----------------------------------[Locals]-----------------------------------!
  integer :: phi_new
  integer :: error
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Global)
!==============================================================================!

  call Mpi_Allreduce(phi,             & ! send buffer
                     phi_new,         & ! recv buffer
                     1,               & ! length
                     comm_type_int,   & ! datatype
                     MPI_SUM,         & ! operation
                     MPI_COMM_WORLD,  &
                     error)

  phi = phi_new

  end subroutine
