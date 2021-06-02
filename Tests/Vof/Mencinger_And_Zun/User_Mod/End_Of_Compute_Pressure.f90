!==============================================================================!
  subroutine User_Mod_End_Of_Compute_Pressure(Flow, Vof, Sol, curr_dt, ini)
!------------------------------------------------------------------------------!
!   This function is called at the end of Compute_Pressure function.           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),  target :: Flow
  type(Vof_Type),    target :: Vof
  type(Solver_Type), target :: Sol
  integer, intent(in)       :: curr_dt  ! current time step
  integer, intent(in)       :: ini      ! inner iteration
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),   pointer :: Grid
  type(Var_Type),    pointer :: p, pp
  type(Matrix_Type), pointer :: a
  real, contiguous,  pointer :: b(:)
  integer                    :: c
  character(22)              :: file_name = 'pressure_xxxxx_yyy.dat'
!==============================================================================!

  ! Take aliases
  Grid => Flow % pnt_grid
  p    => Flow % p
  pp   => Flow % pp
  a    => Sol  % a
  b    => Sol  % b % val

  write(file_name(10:14), '(i5.5)') curr_dt
  write(file_name(16:18), '(i3.3)') ini

  open(99, file=file_name)
  write(99, '(a)') '# compute_pressure: 1:x, 2:p, 3:pp, 4:p%x, 5:pp%x ' //  &
                   ' 4:density  5:diag  6:b'

  do c = 1, Grid % n_cells
    if(Math % Approx_Real(Grid % yc(c), 0.0) .and.  &
       Math % Approx_Real(Grid % zc(c), 0.0)) then
      write(99, '(99es15.5)')  &
        Grid % xc(c), p % n(c), pp % n(c), p % x(c), pp % x(c),  &
        Flow % density(c), a % val(a % dia(c)), b(c)
    end if
  end do

  close(99)

  end subroutine
