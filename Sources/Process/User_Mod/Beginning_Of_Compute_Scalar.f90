!==============================================================================!
  subroutine User_Mod_Beginning_Of_Compute_Scalar(flow, turb, mult, sol,  &
                                                  curr_dt, ini, sc)
!------------------------------------------------------------------------------!
!   This function is called at the end of Compute_Scalar function.             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),      target :: flow
  type(Turb_Type),       target :: turb
  type(Multiphase_Type), target :: mult
  type(Solver_Type),     target :: sol
  integer, intent(in)           :: curr_dt  ! current time step
  integer, intent(in)           :: ini      ! inner iteration
  integer, intent(in)           :: sc       ! scalar index
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),   pointer :: grid
  type(Var_Type),    pointer :: phi, p
  type(Matrix_Type), pointer :: m
!==============================================================================!

  ! Take aliases
  grid => flow % pnt_grid
  p    => flow % p
  phi  => flow % scalar(sc)
  m    => sol % m

  end subroutine
