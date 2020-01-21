include '../User_Mod/Impinging_Jet_Nu.f90'
include '../User_Mod/Impinging_Jet_Profiles.f90'

!==============================================================================!
  subroutine User_Mod_Save_Results(flow, turb, mult, swarm, ts)
!------------------------------------------------------------------------------!
!   Calls User_Impinging_Jet_Nu and User_Impinging_Jet_Profile functions.      !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type)      :: flow
  type(Turb_Type)       :: turb
  type(Multiphase_Type) :: mult
  type(Swarm_Type)      :: swarm
  integer               :: ts     ! time step
!==============================================================================!

  call User_Mod_Impinging_Jet_Nu      (turb)
  call User_Mod_Impinging_Jet_Profiles(turb)

  end subroutine
