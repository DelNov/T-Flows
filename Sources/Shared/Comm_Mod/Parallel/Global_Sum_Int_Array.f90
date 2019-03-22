!==============================================================================!
  subroutine Comm_Mod_Global_Sum_Int_Array(n, phi)
!------------------------------------------------------------------------------!
!   Estimates global sum over all processors.                                  !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer :: n
  integer :: phi(n)
!-----------------------------------[Locals]-----------------------------------!
  integer, allocatable :: phi_res(:)
  integer              :: error
!==============================================================================!

  allocate(phi_res(n))

  call Mpi_Allreduce(phi,             & ! send buffer
                     phi_res,         & ! recv buffer
                     n,               & ! length
                     MPI_INTEGER8,    & ! datatype
                     MPI_SUM,         & ! operation
                     MPI_COMM_WORLD,  &
                     error)

  phi(1:n) = phi_res(1:n)

  end subroutine
