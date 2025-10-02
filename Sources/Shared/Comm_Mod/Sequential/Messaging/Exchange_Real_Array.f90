!==============================================================================!
  subroutine Exchange_Real_Array(Comm, length, phi, dest)
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
  class(Comm_Type), intent(in)    :: Comm         !! communication object
  integer,          intent(in)    :: length       !! length of the array
  real,             intent(inout) :: phi(length)  !! array to be exchanged
  integer,          intent(in)    :: dest         !! destination processor
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Comm)
  Unused(length)
  Unused(phi)
  Unused(dest)
!==============================================================================!

  end subroutine
