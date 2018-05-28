!==============================================================================!
  subroutine Comm_Mod_Global_Max_Real(phi) 
!------------------------------------------------------------------------------!
!   Estimates global maximum among all processors.                             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real    :: phi
!-----------------------------------[Locals]-----------------------------------!
  real    :: phi_new
  integer :: error
!==============================================================================!

  call Mpi_Allreduce(phi,                   & ! send buffer
                     phi_new,               & ! recv buffer 
                     1,                     & ! length     
                     MPI_DOUBLE_PRECISION,  & ! datatype  
                     MPI_MAX,               & ! operation 
                     MPI_COMM_WORLD,        &             
                     error) 

  phi = phi_new

  end subroutine
