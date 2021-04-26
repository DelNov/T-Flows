!==============================================================================!
  subroutine User_Mod_End_Of_Compute_Pressure(flow, mult, curr_dt, ini)
!------------------------------------------------------------------------------!
!   This function is called at the end of Compute_Pressure function.           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),      target :: flow
  type(Multiphase_Type), target :: mult
  integer, intent(in)           :: curr_dt  ! current time step
  integer, intent(in)           :: ini      ! inner iteration
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: grid
  type(Var_Type),  pointer :: p, pp
  integer                  :: c
  character(30)            :: file_name = 'compute_pressure_xxxxx_yyy.dat'
!==============================================================================!

  ! Take aliases
  grid => flow % pnt_grid
  p    => flow % p
  pp   => flow % pp

  write(file_name(18:22), '(i5.5)') curr_dt
  write(file_name(24:26), '(i3.3)') ini

  open(99, file=file_name)
  write(99, '(a)') '# compute_pressure: 1:x, 2:rho, 3:pp, 4:p'

  do c = 1, grid % n_cells
    if(Math_Mod_Approx_Real(grid % yc(c), 0.0) .and.  &
       Math_Mod_Approx_Real(grid % zc(c), 0.0)) then
      write(99, '(99e15.5)')  &
        grid % xc(c), pp % n(c), p % n(c), flow % density (c)
    end if
  end do

  close(99)

  end subroutine
