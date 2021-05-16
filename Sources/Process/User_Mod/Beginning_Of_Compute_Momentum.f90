!==============================================================================!
  subroutine User_Mod_Beginning_Of_Compute_Momentum(Flow, turb, Vof, Sol,  &
                                                    curr_dt, ini)
!------------------------------------------------------------------------------!
!   This function is called at the beginning of Compute_Momentum function.     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),    target :: Flow
  type(Turb_Type),     target :: turb
  type(Vof_Type),      target :: Vof
  type(Solver_Type),   target :: Sol
  integer, intent(in)         :: curr_dt  ! current time step
  integer, intent(in)         :: ini      ! inner iteration
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),   pointer :: grid
  type(Var_Type),    pointer :: u, v, w, p
  type(Matrix_Type), pointer :: M
!==============================================================================!

  ! Take aliases
  grid => Flow % pnt_grid
  u    => Flow % u
  v    => Flow % u
  w    => Flow % u
  p    => Flow % p
  M    => Sol % M

  end subroutine
