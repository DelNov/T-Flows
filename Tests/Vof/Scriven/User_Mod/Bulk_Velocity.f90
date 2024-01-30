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
  type(Bulk_Type), pointer :: bulk
  real                     :: velocity_out
!==============================================================================!

  ! Take alias
  bulk => Flow % bulk

  ! Set some initial value
  bulk_vel = 0.0

  ! By default, assume all outflow is defined as convective
  if(Flow % mass_transfer_model .ne. NO_MASS_TRANSFER) then

    ! vol_src: volume source from mass transfer
    velocity_out = bulk % vol_src / bulk % area_out
  else

    ! vol_in: volume source from inlet velocity
    velocity_out = bulk % vol_in / bulk % area_out
  end if

  bulk_vel = velocity_out

  end subroutine
