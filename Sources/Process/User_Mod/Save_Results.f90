!==============================================================================!
  subroutine User_Mod_Save_Results(flow, turb, mult, swarm, time_step)
!------------------------------------------------------------------------------!
!   This subroutine is called each                                             !
!     RESULTS_SAVE_INTERVAL (set in conrtol file),                             !
!     at the end of a simulation                                               !
!     and after 'save_now' command                                             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),      target :: flow
  type(Turb_Type),       target :: turb
  type(Multiphase_Type), target :: mult
  type(Swarm_Type),      target :: swarm
  integer                       :: time_step
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: grid
  type(Bulk_Type), pointer :: bulk
  type(Var_Type),  pointer :: u, v, w, t
!==============================================================================!

  ! Take aliases
  grid => flow % pnt_grid
  bulk => flow % bulk
  u    => flow % u
  v    => flow % v
  w    => flow % w
  t    => flow % t

  end subroutine
