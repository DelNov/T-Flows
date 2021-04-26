!==============================================================================!
  subroutine User_Mod_End_Of_Compute_Momentum(flow, turb, mult, curr_dt, ini)
!------------------------------------------------------------------------------!
!   This function is called at the end of Compute_Momentum function.           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),      target :: flow
  type(Turb_Type),       target :: turb
  type(Multiphase_Type), target :: mult
  integer, intent(in)           :: curr_dt  ! current time step
  integer, intent(in)           :: ini      ! inner iteration
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: grid
  type(Var_Type),  pointer :: u, p
  integer                  :: c
  character(30)            :: file_name = 'compute_momentum_xxxxx_yyy.dat'
!==============================================================================!

  ! Take aliases
  grid => flow % pnt_grid
  u    => flow % u
  p    => flow % p

  write(file_name(18:22), '(i5.5)') curr_dt
  write(file_name(24:26), '(i3.3)') ini

  open(99, file=file_name)
  write(99, '(a)') '#  User_Mod_End_Of_Compute_Momentum '         //  &
                   ' 1:x  2:u  3:rho, 4:dp/dx*dv'
  do c = 1, grid % n_cells
    if(Math_Mod_Approx_Real(grid % yc(c), 0.0) .and.  &
       Math_Mod_Approx_Real(grid % zc(c), 0.0)) then
      write(99, '(99e15.5)')  &
        grid % xc(c), u % n(c), flow % density(c), p % x(c) * grid % vol(c)
    end if
  end do

  close(99)

  end subroutine
