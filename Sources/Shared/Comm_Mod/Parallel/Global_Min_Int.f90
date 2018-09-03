!==============================================================================!
  subroutine Comm_Mod_Global_Min_Int(phi)
!------------------------------------------------------------------------------!
!   Estimates global min among all processors.                                 !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer :: phi
!-----------------------------------[Locals]-----------------------------------!
  integer :: phi_new
  integer :: error
!==============================================================================!

  call Mpi_Allreduce(phi,             & ! send buffer
                     phi_new,         & ! recv buffer
                     1,               & ! length
                     MPI_INTEGER8,    & ! datatype
                     MPI_MIN,         & ! operation
                     MPI_COMM_WORLD,  &
                     error)

  phi = phi_new

  end subroutine
