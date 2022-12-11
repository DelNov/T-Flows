!==============================================================================!
  subroutine Comm_Mod_Global_Lor_Log(phi)
!------------------------------------------------------------------------------!
!   Estimates logical or over all processors.                                  !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  logical :: phi
!-----------------------------------[Locals]-----------------------------------!
  logical :: phi_new
  integer :: error
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
