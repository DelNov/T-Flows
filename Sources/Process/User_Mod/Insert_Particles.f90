!==============================================================================!
  subroutine User_Mod_Insert_Particles(flow, turb, mult, swarm, n, time)
!------------------------------------------------------------------------------!
!   Insert particles for Lagrangian particle tracking                          !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),      target :: flow
  type(Turb_Type),       target :: turb
  type(Multiphase_Type), target :: mult
  type(Swarm_Type),      target :: swarm
  integer                       :: n     ! time step
  real                          :: time  ! physical time
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: grid
  type(Var_Type),  pointer :: u, v, w, t, phi, vof
!==============================================================================!

  ! Take aliases
  grid => flow % pnt_grid

  end subroutine
