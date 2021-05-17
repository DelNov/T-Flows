!==============================================================================!
  subroutine User_Mod_End_Of_Compute_Momentum(Flow, turb, Vof, Sol,  &
                                               curr_dt, ini)
!------------------------------------------------------------------------------!
!   This function is called at the end of Compute_Momentum function.           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),  target :: Flow
  type(Turb_Type),   target :: turb
  type(Vof_Type),    target :: Vof
  type(Solver_Type), target :: Sol
  integer, intent(in)       :: curr_dt  ! current time step
  integer, intent(in)       :: ini      ! inner iteration
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),   pointer :: grid
  type(Var_Type),    pointer :: u, p
  type(Matrix_Type), pointer :: m
  integer                    :: c
  character(40)              :: file_name =  &
                                'velocity_before_correction_xxxxx_yyy.dat'
!==============================================================================!

  ! Take aliases
  grid => Flow % pnt_grid
  u    => Flow % u
  p    => Flow % p
  m    => Sol  % m

  write(file_name(28:32), '(i5.5)') curr_dt
  write(file_name(34:36), '(i3.3)') ini

  open(99, file=file_name)
  write(99, '(a)') '#  User_Mod_End_Of_Compute_Momentum '         //  &
                   ' 1:x  2:u  3:rho, 4:dp/dx*dv,  5:diag(matrix),  6:b'
  do c = 1, grid % n_cells
    if(Math_Mod_Approx_Real(grid % yc(c), 0.0) .and.  &
       Math_Mod_Approx_Real(grid % zc(c), 0.0)) then
      write(99, '(99es15.5)')  &
        grid % xc(c), u % n(c), Flow % density(c), p % x(c) * grid % vol(c),  &
        m % val(m % dia(c)), Flow % fx(c) - p % x(c) * grid % vol(c),  &
        grid % vol(c)
    end if
  end do

  close(99)

  end subroutine
