!==============================================================================!
  subroutine User_Mod_End_Of_Compute_Momentum(Flow, Turb, Vof, Sol)
!------------------------------------------------------------------------------!
!   This function is called at the end of Compute_Momentum function.           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),  target :: Flow
  type(Turb_Type),   target :: Turb
  type(Vof_Type),    target :: Vof
  type(Solver_Type), target :: Sol
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),   pointer :: Grid
  type(Var_Type),    pointer :: u, p
  type(Matrix_Type), pointer :: M
  integer                    :: c
  character(40)              :: file_name =  &
                                'velocity_before_correction_xxxxx_yyy.dat'
!==============================================================================!

  ! Take aliases
  Grid => Flow % pnt_grid
  u    => Flow % u
  p    => Flow % p
  M    => Sol % Nat % M

  write(file_name(28:32), '(i5.5)') Time % Curr_Dt()
  write(file_name(34:36), '(i3.3)') Iter % Current()

  open(99, file=file_name)
  write(99, '(a)') '#  User_Mod_End_Of_Compute_Momentum '         //  &
                   ' 1:x  2:u  3:rho, 4:dp/dx*dv,  5:diag(matrix),  6:b'
  do c = Cells_In_Domain_And_Buffers()
    if(Math % Approx_Real(Grid % yc(c), 0.0) .and.  &
       Math % Approx_Real(Grid % zc(c), 0.0)) then
      write(99, '(99es15.5)')  &
        Grid % xc(c), u % n(c), Flow % density(c), p % x(c) * Grid % vol(c),  &
        M % val(M % dia(c)), Flow % fx(c) - p % x(c) * Grid % vol(c),  &
        Grid % vol(c)
    end if
  end do

  close(99)

  end subroutine
