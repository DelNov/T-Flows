!==============================================================================!
  subroutine Comm_Mod_Global_Lor_Log_Array(n, phi)
!------------------------------------------------------------------------------!
!   Estimates logical or over all processors.                                  !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer :: n
  logical :: phi(n)
!-----------------------------------[Locals]-----------------------------------!
  logical, allocatable :: phi_res(:)
  integer              :: error
!==============================================================================!

  allocate(phi_res(n))

  call Mpi_Allreduce(phi,             & ! send buffer
                     phi_res,         & ! recv buffer
                     n,               & ! length
                     MPI_LOGICAL8,    & ! datatype
                     MPI_LOR,         & ! operation
                     MPI_COMM_WORLD,  &
                     error)

  phi(1:n) = phi_res(1:n)

  end subroutine
