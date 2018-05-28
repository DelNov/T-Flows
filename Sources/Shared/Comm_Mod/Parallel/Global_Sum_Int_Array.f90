!==============================================================================!
  subroutine Comm_Mod_Global_Sum_Int_Array(phi, n) 
!------------------------------------------------------------------------------!
!   Estimates global summ among all processors.                                !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer :: n
  integer :: phi(n)
!-----------------------------------[Locals]-----------------------------------!
  integer :: phi_new(n)
  integer :: error
!==============================================================================!

  call Mpi_Allreduce(phi,             & ! send buffer
                     phi_new,         & ! recv buffer 
                     n,               & ! length     
                     MPI_INTEGER8,    & ! datatype  
                     MPI_SUM,         & ! operation 
                     MPI_COMM_WORLD,  &             
                     error) 

  phi = phi_new

  end subroutine
