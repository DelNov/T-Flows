!==============================================================================!
  subroutine User_Mod_End_Of_Compute_Pressure(flow, mult, sol, curr_dt, ini)
!------------------------------------------------------------------------------!
!   This function is called at the end of Compute_Pressure function.           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),      target :: flow
  type(Multiphase_Type), target :: mult
  type(Solver_Type),     target :: sol
  integer, intent(in)           :: curr_dt  ! current time step
  integer, intent(in)           :: ini      ! inner iteration
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),   pointer :: grid
  type(Var_Type),    pointer :: u, v, w, p
  type(Matrix_Type), pointer :: a, m        ! pressure and momentum matrices
!==============================================================================!

  ! Take aliases
  grid => flow % pnt_grid
  u    => flow % u
  v    => flow % u
  w    => flow % u
  p    => flow % p
  a    => sol % a
  m    => sol % m

  end subroutine
