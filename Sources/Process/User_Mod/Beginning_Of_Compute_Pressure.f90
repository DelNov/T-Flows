!==============================================================================!
  subroutine User_Mod_Beginning_Of_Compute_Pressure(flow, mult,  &
                                                    curr_dt, ini)
!------------------------------------------------------------------------------!
!   This function is called at the beginning of Compute_Pressure function.     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),      target :: flow
  type(Multiphase_Type), target :: mult
  integer, intent(in)           :: curr_dt  ! current time step
  integer, intent(in)           :: ini      ! inner iteration
!==============================================================================!

  end subroutine
