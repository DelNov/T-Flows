!==============================================================================!
  subroutine User_Mod_Beginning_Of_Compute_Vof(Vof, Sol)
!------------------------------------------------------------------------------!
!   This function is called at the beginning of Compute_Vof function.          !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Vof_Type),    target :: Vof
  type(Solver_Type), target :: Sol
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
  A    => Sol % Nat % A

  if(Time % Curr_Dt() > 120) then
    ! Do something

  end if

  end subroutine
