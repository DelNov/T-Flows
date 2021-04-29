!==============================================================================!
  subroutine User_Mod_Beginning_Of_Compute_Energy(flow, turb, mult, sol,  &
                                                  curr_dt, ini)
!------------------------------------------------------------------------------!
!   This function is called at the beginning of Compute_Energy function.       !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),      target :: flow
  type(Turb_Type),       target :: turb
  type(Multiphase_Type), target :: mult
  type(Solver_Type),     target :: sol
  integer, intent(in)           :: curr_dt  ! current time step
  integer, intent(in)           :: ini      ! inner iteration
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),   pointer :: grid
  type(Var_Type),    pointer :: t, p
  type(Matrix_Type), pointer :: a
!==============================================================================!

  ! Take aliases
  grid => flow % pnt_grid
  t    => flow % t
  p    => flow % p
  a    => sol % a

  end subroutine
