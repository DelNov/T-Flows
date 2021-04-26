!==============================================================================!
  subroutine User_Mod_Beginning_Of_Compute_Scalar(flow, turb, mult,  &
                                                  curr_dt, ini)
!------------------------------------------------------------------------------!
!   This function is called at the end of Compute_Scalar function.             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),      target :: flow
  type(Turb_Type),       target :: turb
  type(Multiphase_Type), target :: mult
  integer, intent(in)           :: curr_dt  ! current time step
  integer, intent(in)           :: ini      ! inner iteration
!==============================================================================!

  end subroutine
