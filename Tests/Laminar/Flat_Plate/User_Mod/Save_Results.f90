# ifdef __INTEL_COMPILER
#   include "User_Mod/Plain_Profiles.f90"
#   include "User_Mod/Plain_Nu.f90"
# else
#   include "Plain_Profiles.f90"
#   include "Plain_Nu.f90"
# endif

!==============================================================================!
  subroutine User_Mod_Save_Results(Flow, Turb, Vof, Swarm, domain)
!------------------------------------------------------------------------------!
!   Calls user-define subroutines                                              !
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

  call User_Mod_Plain_Profiles(Flow, Turb)
  call User_Mod_Plain_Nu      (Flow, Turb)

  end subroutine
