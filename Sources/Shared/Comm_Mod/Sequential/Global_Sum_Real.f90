!==============================================================================!
  subroutine Sum_Real(Global, phi)
!------------------------------------------------------------------------------!
!   Dummy function for sequential runs.                                        !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type), intent(in)    :: Global
  real,             intent(inout) :: phi
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Global)
  Unused(phi)
!==============================================================================!

  end subroutine
