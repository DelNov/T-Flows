!==============================================================================!
  subroutine User_Mod_Beginning_Of_Compute_Vof(Vof, Sol, curr_dt)
!------------------------------------------------------------------------------!
!   This function is called at the beginning of Compute_Vof function.          !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Vof_Type),      target :: Vof
  type(Solver_Type),   target :: Sol
  integer, intent(in)         :: curr_dt  ! current time step
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),   pointer :: Grid
  type(Field_Type),  pointer :: Flow
  type(Var_Type),    pointer :: fun
  type(Matrix_Type), pointer :: A
!==============================================================================!

  ! Take aliases
  Flow => Vof  % pnt_flow
  Grid => Flow % pnt_grid
  fun  => Vof % fun

  end subroutine
