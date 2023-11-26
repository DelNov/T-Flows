!==============================================================================!
  subroutine Lor_Log(Global, phi)
!------------------------------------------------------------------------------!
!>  Performs a logical OR operation across all processors for a logical value.
!>  This subroutine is used in parallel computing within the MPI framework to
!>  ensure all processors reach a common logical state.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type), intent(in)    :: Global  !! global communication object
  logical,          intent(inout) :: phi     !! operand
!-----------------------------------[Locals]-----------------------------------!
  logical :: phi_new
  integer :: error
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Global)
!==============================================================================!

  call Mpi_Allreduce(phi,             & ! send buffer
                     phi_new,         & ! recv buffer
                     1,               & ! length
                     comm_type_log,   & ! datatype
                     MPI_LOR,         & ! operation
                     MPI_COMM_WORLD,  &
                     error)

  phi = phi_new

  end subroutine
