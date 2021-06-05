!==============================================================================!
  subroutine User_Mod_End_Of_Compute_Momentum(Flow, turb, Vof, Sol,  &
                                              curr_dt, ini)
!------------------------------------------------------------------------------!
!   This function is called at the end of Compute_Momentum function.           !
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
  type(Grid_Type),   pointer :: Grid
  type(Var_Type),    pointer :: u, v, w, p
  type(Matrix_Type), pointer :: M
  integer                    :: c
!==============================================================================!

  ! Take aliases
  Grid => Flow % pnt_grid
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
