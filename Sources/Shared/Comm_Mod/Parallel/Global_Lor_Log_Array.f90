!==============================================================================!
  subroutine Comm_Mod_Global_Lor_Log_Array(n, phi)
!------------------------------------------------------------------------------!
!   Estimates logical or over all processors.                                  !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer, intent(in)    :: n
  logical, intent(inout) :: phi(n)
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

  end subroutine
