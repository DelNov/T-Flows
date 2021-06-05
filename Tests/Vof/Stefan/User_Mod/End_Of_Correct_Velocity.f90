!==============================================================================!
  subroutine User_Mod_End_Of_Correct_Velocity(Flow, Vof, Sol, curr_dt, ini)
!------------------------------------------------------------------------------!
!   This function is called at the end of Correct_Velocity function.           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),    target :: Flow
  type(Vof_Type),      target :: Vof
  type(Solver_Type),   target :: Sol
  integer, intent(in)         :: curr_dt
  integer, intent(in)         :: ini
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),   pointer :: grid
  type(Var_Type),    pointer :: u, v, w, p
  type(Matrix_Type), pointer :: M
  integer                    :: c
!==============================================================================!

  ! Take aliases
  grid => Flow % pnt_grid
  u    => Flow % u
  v    => Flow % u
  w    => Flow % u
  p    => Flow % p
  M    => Sol % M

!@  ! Nullify v and w velocity components
!@  do c = -Grid % n_bnd_cells, Grid % n_cells
!@    v % n(c) = 0
!@    w % n(c) = 0
!@  end do

  end subroutine
