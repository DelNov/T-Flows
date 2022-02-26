!==============================================================================!
  subroutine Alias_T2(Turb, t2)
!------------------------------------------------------------------------------!
!   Create alias for t2 variable.                                              !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Turb_Type), target  :: Turb
  type(Var_Type),   pointer :: t2
!==============================================================================!

  t2 => Turb % t2

  end subroutine
