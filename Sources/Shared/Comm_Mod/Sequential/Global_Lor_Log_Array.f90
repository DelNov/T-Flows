!==============================================================================!
  subroutine Lor_Log_Array(Global, n, phi)
!------------------------------------------------------------------------------!
!   Dummy function for sequential runs.                                        !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type), intent(in)    :: Global
  integer,          intent(in)    :: n
  logical,          intent(inout) :: phi(n)
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Global)
  Unused(n)
  Unused(phi)
!==============================================================================!

  end subroutine
