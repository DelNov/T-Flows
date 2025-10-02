!==============================================================================!
  subroutine User_Mod_End_Of_Compute_Momentum(Flow, Turb, Vof, Sol)
!------------------------------------------------------------------------------!
!   This function is called at the end of Compute_Momentum function.           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),  target :: Flow
  type(Turb_Type),   target :: Turb
  type(Vof_Type),    target :: Vof
  type(Solver_Type), target :: Sol
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),   pointer :: Grid
  type(Var_Type),    pointer :: u, v, w, p
  type(Matrix_Type), pointer :: M
!==============================================================================!

  ! Take aliases
  Grid => Flow % pnt_grid
  u    => Flow % u
  v    => Flow % u
  w    => Flow % u
  p    => Flow % p
  M    => Sol % Nat % M

  end subroutine
