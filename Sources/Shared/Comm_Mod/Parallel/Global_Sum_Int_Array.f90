!==============================================================================!
  subroutine Comm_Mod_Global_Sum_Int_Array(n, phi)
!------------------------------------------------------------------------------!
!   Estimates global sum over all processors.                                  !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer, intent(in)    :: n
  integer, intent(inout) :: phi(n)
!-----------------------------------[Locals]-----------------------------------!
  integer :: error
!==============================================================================!

  call Mpi_Allreduce(MPI_IN_PLACE,    & ! indicate that send and recv are same
                     phi,             & ! send and recv buffer
                     n,               & ! length
                     comm_type_int,   & ! datatype
                     MPI_SUM,         & ! operation
                     MPI_COMM_WORLD,  &
                     error)

  end subroutine
