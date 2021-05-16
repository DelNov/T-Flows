!==============================================================================!
  subroutine User_Mod_Save_Results(Flow, turb, Vof, swarm, ts)
!------------------------------------------------------------------------------!
!   This subroutine is called each RESULTS_SAVE_INTERVAL (set in control       !
!   file), at the end of a simulation and after 'save_now' command.            !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),    target :: Flow
  type(Turb_Type),     target :: turb
  type(Vof_Type),      target :: Vof
  type(Swarm_Type),    target :: swarm
  integer, intent(in)         :: ts     ! time step
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: grid
  type(Bulk_Type), pointer :: bulk
  type(Var_Type),  pointer :: u, v, w, t
!==============================================================================!

  ! Take aliases
  grid => Flow % pnt_grid
  bulk => Flow % bulk
  u    => Flow % u
  v    => Flow % v
  w    => Flow % w
  t    => Flow % t

  end subroutine
