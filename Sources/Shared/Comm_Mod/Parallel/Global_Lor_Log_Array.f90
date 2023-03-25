!==============================================================================!
  subroutine Lor_Log_Array(Global, n, phi)
!------------------------------------------------------------------------------!
!   Estimates logical or over all processors.                                  !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type), intent(in)    :: Global
  integer,          intent(in)    :: n
  logical,          intent(inout) :: phi(n)
!-----------------------------------[Locals]-----------------------------------!
  integer :: error
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Global)
!==============================================================================!

  call Mpi_Allreduce(MPI_IN_PLACE,    & ! indicate that send and recv are same
                     phi,             & ! send and recv buffer
                     n,               & ! length
                     comm_type_log,   & ! datatype
                     MPI_LOR,         & ! operation
                     MPI_COMM_WORLD,  &
                     error)

  end subroutine
