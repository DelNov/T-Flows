!==============================================================================!
  subroutine Max_Int(Global, phi)
!------------------------------------------------------------------------------!
!>  The subroutine searches for maximum integer value across all processors.
!>  This subroutine is useful when a global maximum needs to be determined from
!>  values distributed across multiple processors in a parallel computation.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type), intent(in)    :: Global  !! global communication object
  integer,          intent(inout) :: phi     !! maximum over all processors
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
