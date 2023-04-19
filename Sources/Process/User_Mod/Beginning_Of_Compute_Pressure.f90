!==============================================================================!
  subroutine User_Mod_Beginning_Of_Compute_Pressure(Flow, Vof, Sol)
!------------------------------------------------------------------------------!
!   This function is called at the beginning of Compute_Pressure function.     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),  target :: Flow
  type(Vof_Type),    target :: Vof
  type(Solver_Type), target :: Sol
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),   pointer :: Grid
  type(Var_Type),    pointer :: u, v, w, p
  type(Matrix_Type), pointer :: A, M        ! pressure and momentum matrices
!==============================================================================!

  ! Take aliases
  Grid => Flow % pnt_grid
  u    => Flow % u
  v    => Flow % u
  w    => Flow % u
  p    => Flow % p
  A    => Sol % Nat % A
  M    => Sol % Nat % M

  end subroutine
