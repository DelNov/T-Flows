!==============================================================================!
  subroutine User_Mod_End_Of_Compute_Energy(flow, turb, mult, ini)
!------------------------------------------------------------------------------!
!   This function is called at the end of Compute_Energy function.             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),      target :: flow
  type(Turb_Type),       target :: turb
  type(Multiphase_Type), target :: mult
  integer, intent(in)           :: ini   ! inner iteration
!==============================================================================!

  end subroutine
