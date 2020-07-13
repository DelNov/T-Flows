!==============================================================================!
  subroutine User_Mod_End_Of_Correct_Velocity(flow, mult, sol, dt, ini)
!------------------------------------------------------------------------------!
!   This function is called at the end of Correct_Velocity function.           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),      target :: flow
  type(Multiphase_Type), target :: mult
  type(Solver_Type),     target :: sol
  real                          :: dt
  integer                       :: ini
!==============================================================================!

  end subroutine
