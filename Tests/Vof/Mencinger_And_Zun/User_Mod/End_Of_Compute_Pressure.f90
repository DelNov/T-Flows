!==============================================================================!
  subroutine User_Mod_End_Of_Compute_Pressure(Flow, Vof, Sol)
!------------------------------------------------------------------------------!
!   This function is called at the end of Compute_Pressure function.           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),  target :: Flow
  type(Vof_Type),    target :: Vof
  type(Solver_Type), target :: Sol
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),   pointer :: Grid
  type(Var_Type),    pointer :: p, pp
  type(Matrix_Type), pointer :: A
  real, contiguous,  pointer :: b(:)
  integer                    :: c
  character(22)              :: file_name = 'pressure_xxxxx_yyy.dat'
!==============================================================================!

  ! Take aliases
  Grid => Flow % pnt_grid
  p    => Flow % p
  pp   => Flow % pp
  A    => Sol % Nat % A
  b    => Sol % Nat % b % val

  write(file_name(10:14), '(i5.5)') Time % Curr_Dt()
  write(file_name(16:18), '(i3.3)') Iter % Current()

  open(99, file=file_name)
  write(99, '(A)') '# compute_pressure: 1:x, 2:p, 3:pp, 4:p%x, 5:pp%x ' //  &
                   ' 4:density  5:diag  6:b'

  do c = Cells_In_Domain_And_Buffers()
    if(Math % Approx_Real(Grid % yc(c), 0.0) .and.  &
       Math % Approx_Real(Grid % zc(c), 0.0)) then
      write(99, '(99es15.5)')  &
        Grid % xc(c), p % n(c), pp % n(c), p % x(c), pp % x(c),  &
        Flow % density(c), A % val(A % dia(c)), b(c)
    end if
  end do

  close(99)

  end subroutine
