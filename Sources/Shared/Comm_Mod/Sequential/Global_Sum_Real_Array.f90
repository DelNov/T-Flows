!==============================================================================!
  subroutine Sum_Real_Array(Global, n, phi)
!------------------------------------------------------------------------------!
!   Dummy function for sequential runs.                                        !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type), intent(in)    :: Global
  integer,          intent(in)    :: n
  real,             intent(inout) :: phi(n)
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Global)
  Unused(n)
  Unused(phi)
!==============================================================================!

  end subroutine
