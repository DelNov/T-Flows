# ifdef __INTEL_COMPILER
#   include "User_Mod/Backstep_Profiles.f90"
#   include "User_Mod/Backstep_Cf_St.f90"
# else
#   include "Backstep_Profiles.f90"
#   include "Backstep_Cf_St.f90"
# endif

!==============================================================================!
  subroutine User_Mod_Save_Results(Flow, Turb, Vof, Swarm, domain)
!------------------------------------------------------------------------------!
!   Calls User_Backstep_Profiles and User_Backstep_Cf_St functions.            !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type), target :: Flow
  type(Turb_Type),  target :: Turb
  type(Vof_Type),   target :: Vof
  type(Swarm_Type), target :: Swarm
  integer, optional        :: domain
!==============================================================================!

  ! Don't save if this is intial condition, nothing is developed yet
  if(Time % Curr_Dt() .eq. 0) return

  call User_Mod_Backstep_Profiles(Flow, Turb)
  call User_Mod_Backstep_Cf_St   (Flow, Turb)

  end subroutine
