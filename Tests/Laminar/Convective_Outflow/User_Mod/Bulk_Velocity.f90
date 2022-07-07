!==============================================================================!
  subroutine User_Mod_Bulk_Velocity(Flow, bulk_vel)
!------------------------------------------------------------------------------!
!   This function is called after bulk velocity is computed in the function    !
!   Convective_Outflow, so that user can specify custom bulk velocity.  This   !
!   is usefull in non-trivial geometries where monitoring planes cannot be     !
!   defined in an obvious way.                                                 !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type) :: Flow
  real             :: bulk_vel
!==============================================================================!

  bulk_vel = 1.0;

  end subroutine
