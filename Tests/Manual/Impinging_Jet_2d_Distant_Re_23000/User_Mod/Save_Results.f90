# ifdef __INTEL_COMPILER
#   include "User_Mod/Save_Impinging_Jet_Nu.f90"
#   include "User_Mod/Save_Impinging_Jet_Profiles.f90"
# else
#   include "Save_Impinging_Jet_Nu.f90"
#   include "Save_Impinging_Jet_Profiles.f90"
# endif

!==============================================================================!
  subroutine User_Mod_Save_Results(Flow, Turb, Vof, Swarm, domain)
!------------------------------------------------------------------------------!
!   Calls Save_Impinging_Jet_Nu and Save_Impinging_Jet_Profile functions.      !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type)  :: Flow
  type(Turb_Type)   :: Turb
  type(Vof_Type)    :: Vof
  type(Swarm_Type)  :: Swarm
  integer, optional :: domain
!==============================================================================!

  ! Don't save if this is intial condition, nothing is developed yet
  if(Time % Curr_Dt() .eq. 0) return

  call Save_Impinging_Jet_Nu      (Turb)
  call Save_Impinging_Jet_Profiles(Turb)

  end subroutine
