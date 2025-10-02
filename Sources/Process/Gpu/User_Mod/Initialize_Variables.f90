!==============================================================================!
  subroutine User_Mod_Initialize_Variables(Grid, Flow, Turb)
!------------------------------------------------------------------------------!
!   User initialization of dependent variables.                                !
!                                                                              !
!   Note: Be aware that this function is called from CPU only.                 !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type),  target :: Grid
  type(Field_Type), target :: Flow
  type(Turb_Type),  target :: Turb
!-----------------------------------[Locals]-----------------------------------!
!==============================================================================!

  end subroutine
