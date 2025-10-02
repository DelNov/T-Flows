!==============================================================================!
  subroutine Max_Real_Array(Global, n, phi)
!------------------------------------------------------------------------------!
!>  The subroutine Max_Real_Array in a parallel computing framework calculates
!>  the global max of real arrays across all processors. This subroutine is
!>  vital in scenarios where an aggregate max of distributed real values is
!>  required for collective decision-making or further computational processes.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type), intent(in)    :: Global  !! global communication object
  integer,          intent(in)    :: n       !! array size
  real,             intent(inout) :: phi(n)  !! operand array
!-----------------------------------[Locals]-----------------------------------!
  integer :: error
!==============================================================================!

  call Mpi_Allreduce(MPI_IN_PLACE,    & ! indicate that send and recv are same
                     phi,             & ! send and recv buffer
                     n,               & ! length
                     comm_type_real,  & ! datatype
                     MPI_MAX,         & ! operation
                     MPI_COMM_WORLD,  &
                     error)

  ! Although the barrier shouldn't be needed here, in some rare
  ! cases (only some small grids with some number of processors,
  ! with certain combinations of compiler/libraries and OS'), it
  ! was a gurantess that all the processors have the correct data
  ! It was introduced only in conjunction with MPI_IN_PLACE.
  call Global % Wait()

  end subroutine
