!==============================================================================!
  subroutine Lor_Log(Global, phi)
!------------------------------------------------------------------------------!
!   Estimates logical or over all processors.                                  !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type), intent(in)    :: Global
  logical,          intent(inout) :: phi
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
