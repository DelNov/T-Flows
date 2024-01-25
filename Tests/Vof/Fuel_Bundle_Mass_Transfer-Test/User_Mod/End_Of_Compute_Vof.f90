!==============================================================================!
  subroutine User_Mod_End_Of_Compute_Vof(Vof, Sol)
!------------------------------------------------------------------------------!
!   This function is called at the end of Compute_Vof function.                !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Vof_Type),      target :: Vof
  type(Solver_Type),   target :: Sol
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),   pointer :: Grid
  type(Field_Type),  pointer :: Flow
  type(Var_Type),    pointer :: fun
  type(Matrix_Type), pointer :: A
  integer :: c
!==============================================================================!

  ! Take aliases
  Flow => Vof  % pnt_flow
  Grid => Flow % pnt_grid
  fun  => Vof % fun

  do c = 1, Grid % n_cells
    if (Grid % zc(c) > 3.75e-2) then
      Vof % fun % n(c) = 1.0
    endif
  end do

  end subroutine
