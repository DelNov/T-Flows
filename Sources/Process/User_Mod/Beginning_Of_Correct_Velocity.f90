!==============================================================================!
  subroutine User_Mod_Beginning_Of_Correct_Velocity(flow, Vof, Sol,  &
                                                    curr_dt, ini)
!------------------------------------------------------------------------------!
!   This function is called at the beginning of Correct_Velocity function.     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),    target :: flow
  type(Vof_Type),      target :: Vof
  type(Solver_Type),   target :: Sol
  integer, intent(in)         :: curr_dt
  integer, intent(in)         :: ini
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),   pointer :: grid
  type(Var_Type),    pointer :: u, v, w, p
  type(Matrix_Type), pointer :: M
!==============================================================================!

  ! Take aliases
  grid => flow % pnt_grid
  u    => flow % u
  v    => flow % u
  w    => flow % u
  p    => flow % p
  M    => Sol % M

  end subroutine
