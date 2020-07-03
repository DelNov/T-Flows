!==============================================================================!
  subroutine User_Mod_Save_Swarm(flow, turb, mult, swarm, ts)
!------------------------------------------------------------------------------!
!   This subroutine is called each RESULTS_SAVE_INTERVAL (set in control       !
!   file), at the end of a simulation and after 'save_now' command.            !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),      target :: flow
  type(Turb_Type),       target :: turb
  type(Multiphase_Type), target :: mult
  type(Swarm_Type),      target :: swarm
  integer                       :: ts
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
