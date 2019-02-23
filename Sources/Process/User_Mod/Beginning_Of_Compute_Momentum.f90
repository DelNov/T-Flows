!==============================================================================!
  subroutine User_Mod_Beginning_Of_Compute_Momentum(flow, dt, ini)
!------------------------------------------------------------------------------!
!   This function is called at the beginning of Compute_Momentum function.     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type) :: flow
  real             :: dt    ! time step    
  integer          :: ini   ! inner itteration
!==============================================================================!

  end subroutine
