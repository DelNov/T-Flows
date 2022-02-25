!==============================================================================!
  subroutine User_Mod_End_Of_Compute_Energy(Flow, Turb, Vof, Sol, curr_dt, ini)
!------------------------------------------------------------------------------!
!   This function is called at the end of Compute_Energy function.             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),    target :: Flow
  type(Turb_Type),     target :: Turb
  type(Vof_Type),      target :: Vof
  type(Solver_Type),   target :: Sol
  integer, intent(in)         :: curr_dt  ! current time step
  integer, intent(in)         :: ini      ! inner iteration
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),   pointer :: Grid
  type(Var_Type),    pointer :: t, p
  type(Matrix_Type), pointer :: A
!==============================================================================!

  ! Take aliases
  Grid => Flow % pnt_grid
  t    => Flow % t
  p    => Flow % p
  A    => Sol % Nat % A

  end subroutine
