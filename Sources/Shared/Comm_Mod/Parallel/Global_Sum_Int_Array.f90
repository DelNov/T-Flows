!==============================================================================!
  subroutine Comm_Mod_Global_Sum_Int_Array(phi, n)
!------------------------------------------------------------------------------!
!   Estimates global sum over all processors.                                  !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer :: phi(n)
  integer :: n
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
