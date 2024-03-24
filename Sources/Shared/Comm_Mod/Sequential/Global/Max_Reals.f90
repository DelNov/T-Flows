!==============================================================================!
  subroutine Max_Reals(Global, phi_01, phi_02, phi_03, phi_04,  &
                               phi_05, phi_06, phi_07, phi_08)
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
  class(Comm_Type), intent(in)    :: Global  !! global communication object
  real,             intent(inout) :: phi_01  !! global max 1
  real, optional,   intent(inout) :: phi_02  !! global max 2
  real, optional,   intent(inout) :: phi_03  !! global max 3
  real, optional,   intent(inout) :: phi_04  !! global max 4
  real, optional,   intent(inout) :: phi_05  !! global max 5
  real, optional,   intent(inout) :: phi_06  !! global max 6
  real, optional,   intent(inout) :: phi_07  !! global max 7
  real, optional,   intent(inout) :: phi_08  !! global max 8
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Global)
  Unused(phi_01)
  Unused(phi_02)
  Unused(phi_03)
  Unused(phi_04)
  Unused(phi_05)
  Unused(phi_06)
  Unused(phi_07)
  Unused(phi_08)
!==============================================================================!

  end subroutine
