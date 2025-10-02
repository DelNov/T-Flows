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
  type(Field_Type), target :: Flow
  real                     :: bulk_vel
!-----------------------------------[Locals]-----------------------------------!
! type(Grid_Type), pointer :: Grid
  type(Bulk_Type), pointer :: bulk
  real                     :: velocity_in, velocity_out
!==============================================================================!

  ! Take alias
  bulk => Flow % bulk

  ! Set some initial value
  bulk_vel = 0.0

  ! By default, assume all outflow is defined as convective
  if(bulk % area_in > TINY) then
    velocity_in  = bulk % vol_in / bulk % area_in
    velocity_out = bulk % vol_in / bulk % area_out

    bulk_vel = velocity_out
  end if

  end subroutine
