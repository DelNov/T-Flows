!==============================================================================!
  subroutine Min_Real(Global, phi)
!------------------------------------------------------------------------------!
!   Estimates global minimum among all processors.                             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type), intent(in)    :: Global
  real,             intent(inout) :: phi
!-----------------------------------[Locals]-----------------------------------!
  real    :: phi_new
  integer :: error
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Global)
!==============================================================================!

  call Mpi_Allreduce(phi,                   & ! send buffer
                     phi_new,               & ! recv buffer
                     1,                     & ! length
                     comm_type_real,        & ! datatype
                     MPI_MIN,               & ! operation
                     MPI_COMM_WORLD,        &
                     error)

  phi = phi_new

  end subroutine
