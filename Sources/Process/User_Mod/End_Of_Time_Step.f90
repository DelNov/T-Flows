!==============================================================================!
  subroutine User_Mod_End_Of_Time_Step(flow, swarm, n, time)
!------------------------------------------------------------------------------!
!   This function is called at the end of time step.                           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type), target :: flow
  type(Swarm_Type), target :: swarm
  integer                  :: n     ! time step
  real                     :: time  ! physical time
!----------------------------------[Locals]------------------------------------!
  type(Grid_Type), pointer :: grid
  type(Var_Type),  pointer :: u, v, w
!==============================================================================!

  ! Take aliases
  grid => flow % pnt_grid
  u    => flow % u
  v    => flow % v
  w    => flow % w

  end subroutine
