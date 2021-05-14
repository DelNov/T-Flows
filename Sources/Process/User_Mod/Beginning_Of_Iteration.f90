!==============================================================================!
  subroutine User_Mod_Beginning_Of_Iteration(flow, turb, Vof, swarm, n, time)
!------------------------------------------------------------------------------!
!   This function is called at the beginning of time step.                     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),    target :: flow
  type(Turb_Type),     target :: turb
  type(Vof_Type),      target :: Vof
  type(Swarm_Type),    target :: swarm
  integer, intent(in)         :: n     ! time step
  real,    intent(in)         :: time  ! physical time
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: grid
  type(Var_Type),  pointer :: u, v, w, t, phi
!==============================================================================!

  ! Take aliases
  grid => flow % pnt_grid

  end subroutine
