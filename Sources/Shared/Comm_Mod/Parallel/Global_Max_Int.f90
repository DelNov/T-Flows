!==============================================================================!
  subroutine Max_Int(Global, phi)
!------------------------------------------------------------------------------!
!   Estimates global max among all processors.                                 !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type), intent(in)    :: Global
  integer,          intent(inout) :: phi
!-----------------------------------[Locals]-----------------------------------!
  integer :: phi_new
  integer :: error
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Global)
!==============================================================================!

  call Mpi_Allreduce(phi,             & ! send buffer
                     phi_new,         & ! recv buffer
                     1,               & ! length
                     comm_type_int,   & ! datatype
                     MPI_MAX,         & ! operation
                     MPI_COMM_WORLD,  &
                     error)

  phi = phi_new

  end subroutine
