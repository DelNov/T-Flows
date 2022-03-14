!==============================================================================!
  subroutine User_Mod_Beginning_Of_Iteration(Flow, Turb, Vof, Swarm,  &
                                             curr_dt, time)
!------------------------------------------------------------------------------!
!   This function is called at the beginning of iteration.                     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),    target :: Flow
  type(Turb_Type),     target :: Turb
  type(Vof_Type),      target :: Vof
  type(Swarm_Type),    target :: Swarm
  integer, intent(in)         :: curr_dt  ! time step
  real,    intent(in)         :: time     ! physical time
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
  type(Var_Type),  pointer :: u, v, w, t, phi
!==============================================================================!

  ! Take aliases
  Grid => Flow % pnt_grid

  end subroutine
