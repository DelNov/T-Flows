!==============================================================================!
  subroutine Create_New_Types(Comm)
!------------------------------------------------------------------------------!
!>  Dummy function in sequential runs serving as a placeholder for parallel
!>  computing operations when MPI is not available. The functions is designed
!>  to maintain compatibility and uniformity in function calls between parallel
!>  and sequential versions of the software, essentially doing nothing but
!>  ensuring that the code structure remains consistent regardless of the
!>  compilation mode (parallel or sequential).
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type), intent(inout) :: Comm  !! communicator from grid
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Comm)
!==============================================================================!

  end subroutine
