!==============================================================================!
  subroutine Min_Real(Global, phi)
!------------------------------------------------------------------------------!
!>  The subroutine searches for minimum real value across all processors.
!>  This subroutine is useful when a global minimum needs to be determined from
!>  values distributed across multiple processors in a parallel computation.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type), intent(in)    :: Global  !! global communication object
  real,             intent(inout) :: phi     !! minimum over all processors
!-----------------------------------[Locals]-----------------------------------!
  real    :: phi_new
  integer :: error
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Global)
!==============================================================================!

  call Mpi_Allreduce(phi,             & ! send buffer
                     phi_new,         & ! recv buffer
                     1,               & ! length
                     comm_type_real,  & ! datatype
                     MPI_MIN,         & ! operation
                     MPI_COMM_WORLD,  &
                     error)

  phi = phi_new

  end subroutine