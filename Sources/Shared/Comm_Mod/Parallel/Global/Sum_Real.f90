!==============================================================================!
  subroutine Sum_Real(Global, phi)
!------------------------------------------------------------------------------!
!>  The subroutine Sum_Real in a parallel computing framework calculates
!>  the global sum of a real value across all processors. This subroutine is
!>  vital in scenarios where an aggregate sum of a distributed real value is
!>  required for collective decision-making or further computational processes.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type), intent(in)    :: Global  !! global communication object
  real,             intent(inout) :: phi     !! global sum
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
                     MPI_SUM,         & ! operation
                     MPI_COMM_WORLD,  &
                     error)

  phi = phi_new

  end subroutine
