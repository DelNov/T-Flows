!==============================================================================!
  subroutine User_Mod_Save_Swarm(Flow, Turb, Vof, Swarm, domain)
!------------------------------------------------------------------------------!
!   This subroutine is called each RESULTS_SAVE_INTERVAL (set in control       !
!   file), at the end of a simulation and after 'save_now' command.            !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type), target  :: Flow
  type(Turb_Type),  target  :: Turb
  type(Vof_Type),   target  :: Vof
  type(Swarm_Type), target  :: Swarm
  integer,         optional :: domain
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
  type(Bulk_Type), pointer :: bulk
  type(Var_Type),  pointer :: u, v, w, t
!==============================================================================!

  ! Take aliases
  Grid => Flow % pnt_grid
  bulk => Flow % bulk
  u    => Flow % u
  v    => Flow % v
  w    => Flow % w
  t    => Flow % t

  end subroutine
