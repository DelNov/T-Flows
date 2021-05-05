!==============================================================================!
  subroutine User_Mod_Initialize_Variables(flow, turb, mult, swarm, Sol)
!------------------------------------------------------------------------------!
!   User initialization of dependent variables.                                !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),      target :: flow
  type(Turb_Type),       target :: turb
  type(Multiphase_Type), target :: mult
  type(Swarm_Type),      target :: swarm
  type(Solver_Type),     target :: Sol
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: grid
  type(Var_Type),  pointer :: u, v, w, t, phi, vof
!==============================================================================!

  ! Take aliases
  grid => flow % pnt_grid

  end subroutine
