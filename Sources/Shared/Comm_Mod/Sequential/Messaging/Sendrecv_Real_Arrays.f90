!==============================================================================!
  subroutine Sendrecv_Real_Arrays(Comm, len_s, phi_s,  &
                                        len_r, phi_r, dest)
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
  class(Comm_Type), intent(in) :: Comm          !! communication object
  integer,          intent(in) :: len_s         !! send length
  real,             intent(in) :: phi_s(len_s)  !! send buffer
  integer,          intent(in) :: len_r         !! receive length
  real,             intent(in) :: phi_r(len_r)  !! receive buffer
  integer,          intent(in) :: dest          !! destination processor
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Comm)
  Unused(len_s)
  Unused(phi_s)
  Unused(len_r)
  Unused(phi_r)
  Unused(dest)
!==============================================================================!

  end subroutine
