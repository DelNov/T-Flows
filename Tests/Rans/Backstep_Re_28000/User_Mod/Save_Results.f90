include '../User_Mod/Backstep_Profiles.f90'
include '../User_Mod/Backstep_Cf_St.f90'

!==============================================================================!
  subroutine User_Mod_Save_Results(flow, turb, mult, n)
!------------------------------------------------------------------------------!
!   Calls User_Backstep_Profiles and User_Backstep_Cf_St functions.            !
!------------------------------------------------------------------------------!
  use Field_Mod, only: Field_Type
  use Turb_Mod,  only: Turb_Type
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),       target :: flow
  type(Turb_Type),        target :: turb
  type(Multiphase_Type),  target :: mult
  integer                        :: n
!==============================================================================!

  call User_Mod_Backstep_Profiles(flow, turb)
  call User_Mod_Backstep_Cf_St   (flow, turb)

  end subroutine
