!==============================================================================!
  subroutine Lor_Log_Array(Global, n, phi)
!------------------------------------------------------------------------------!
!>  Performs a logical OR operation across all processors for an array of
!>  logical values. This subroutine is used in parallel computing within the
!>  MPI framework to ensure all processors reach a common logical state for
!>  each element of the array.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type), intent(in)    :: Global  !! global communication object
  integer,          intent(in)    :: n       !! array size
  logical,          intent(inout) :: phi(n)  !! operand array
!-----------------------------------[Locals]-----------------------------------!
  integer :: error
!==============================================================================!

  call Mpi_Allreduce(MPI_IN_PLACE,    & ! indicate that send and recv are same
                     phi,             & ! send and recv buffer
                     n,               & ! length
                     comm_type_log,   & ! datatype
                     MPI_LOR,         & ! operation
                     MPI_COMM_WORLD,  &
                     error)

  ! Although the barrier shouldn't be needed here, in some rare
  ! cases (only some small grids with some number of processors,
  ! with certain combinations of compiler/libraries and OS'), it
  ! was a gurantess that all the processors have the correct data
  ! It was introduced only in conjunction with MPI_IN_PLACE.
  call Global % Wait()

  end subroutine
