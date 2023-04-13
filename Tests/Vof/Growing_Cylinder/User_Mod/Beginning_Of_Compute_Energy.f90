!==============================================================================!
  subroutine User_Mod_Beginning_Of_Compute_Energy(Flow, Turb, Vof, Sol,  &
                                                  curr_dt, ini)
!------------------------------------------------------------------------------!
!   This function is called at the beginning of Compute_Energy function.       !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),    target :: Flow
  type(Turb_Type),     target :: Turb
  type(Vof_Type),      target :: Vof
  type(Solver_Type),   target :: Sol
  integer, intent(in)         :: curr_dt  ! current time step
  integer, intent(in)         :: ini      ! inner iteration
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),   pointer :: Grid
  type(Var_Type),    pointer :: t, p
  type(Matrix_Type), pointer :: A
  integer                    :: c
!==============================================================================!

  ! Take aliases
  Grid => Flow % pnt_grid
  t    => Flow % t
  p    => Flow % p
  A    => Sol % Nat % A

  ! If cell was crossed by interface (was liquid; > 0.5) and is now
  ! vapor (< 0.5) set its temperature to saturation
  do c = Cells_In_Domain_And_Buffers()
    if(Vof % fun % o(c) > 0.5 .and. Vof % fun % n(c) < 0.5) then
      t % n(c)  = 100.0
      t % o(c)  = 100.0
      t % oo(c) = 100.0
    end if
  end do

  end subroutine
