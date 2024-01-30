!==============================================================================!
  subroutine User_Mod_End_Of_Correct_Velocity(Flow, Vof, Sol)
!------------------------------------------------------------------------------!
!   This function is called at the end of Correct_Velocity function.           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),    target :: Flow
  type(Vof_Type),      target :: Vof
  type(Solver_Type),   target :: Sol
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),   pointer :: Grid
  type(Var_Type),    pointer :: u, v, w, p
  type(Matrix_Type), pointer :: M
  real :: vel_max,vel_max_before
  real, parameter :: w_inlet = 0.2
  real, parameter :: vlimit = 20.0
  integer :: c
!==============================================================================!

  ! Take aliases
  Grid => Flow % pnt_grid
  u    => Flow % u
  v    => Flow % v
  w    => Flow % w
  p    => Flow % p
  M    => Sol % Nat % M

  do c = -Grid % n_bnd_cells, Grid % n_cells
    if (Grid % zc(c) > 3.65e-2) then
      u % n(c) = 0.0
      v % n(c) = 0.0
      w % n(c) = w_inlet
    end if
  end do

  ! Calculate velocity magnitude for normalization
  vel_max = MICRO
  do c = -Grid % n_bnd_cells, Grid % n_cells
    vel_max = max(vel_max, sqrt(u % n(c)**2 + v % n(c)**2 + w % n(c)**2))
  end do
  call Global % Max_Real(vel_max)
  vel_max_before = vel_max

  do c = -Grid % n_bnd_cells, Grid % n_cells
    u % n(c) = min ( max( u % n(c), -vlimit), vlimit)
    v % n(c) = min ( max( v % n(c), -vlimit), vlimit)
    w % n(c) = min ( max( w % n(c), -vlimit), vlimit)
  end do

  ! Calculate velocity magnitude for normalization
  vel_max = MICRO
  do c = -Grid % n_bnd_cells, Grid % n_cells
    vel_max = max(vel_max, sqrt(u % n(c)**2 + v % n(c)**2 + w % n(c)**2))
  end do
  call Global % Max_Real(vel_max)
  IF (First_proc()) write(*,'(a,2e10.2)') __FILE__//' vel_max = ',  &
                                          vel_max_before, vel_max

  end subroutine
