!==============================================================================!
  subroutine User_Mod_Beginning_Of_Compute_Scalar(Flow, Turb, Vof, Sol, sc)
!------------------------------------------------------------------------------!
!   This function is called at the end of Compute_Scalar function.             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),    target :: Flow
  type(Turb_Type),     target :: Turb
  type(Vof_Type),      target :: Vof
  type(Solver_Type),   target :: Sol
  integer, intent(in)         :: sc       ! scalar index
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),   pointer :: Grid
  type(Var_Type),    pointer :: phi, p
  type(Matrix_Type), pointer :: M
!==============================================================================!

  ! Take aliases
  Grid => Flow % pnt_grid
  p    => Flow % p
  phi  => Flow % scalar(sc)
  M    => Sol % Nat % M

  end subroutine
