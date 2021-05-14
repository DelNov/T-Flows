!==============================================================================!
  subroutine User_Mod_Insert_Particles(flow, turb, Vof, swarm, n, time)
!------------------------------------------------------------------------------!
!   Insert particles for Lagrangian particle tracking                          !
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
