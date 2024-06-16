!==============================================================================!
  subroutine Turb_Update_Host(Gpu, Turb, Flow)
!------------------------------------------------------------------------------!
!>  Copy all the turbulence variables you need for post-processing back to CPU
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Gpu_Type)  :: Gpu   !! parent class
  type(Turb_Type)  :: Turb  !! to check if wall distance should be updated
  type(Field_Type) :: Flow  !! field to transfer to device
!-----------------------[Avoid unused argument warning]------------------------!
# if T_FLOWS_GPU == 0
    Unused(Gpu)
    Unused(Turb)
    Unused(Flow)
# endif
!==============================================================================!

# if T_FLOWS_GPU == 0
    return
# else
    O_Print '(a)', ' # Copying turbulence variables back to the host'
# endif

  !-----------------------!
  !   Smagorinsky model   !
  !-----------------------!
  if(Turb % model .eq. LES_SMAGORINSKY) then
    call Gpu % Vector_Update_Host(Turb % vis_t)
  end if

  !------------------------------------------------!
  !   Variables needed for all turbulence models   !
  !------------------------------------------------!
  if(Turb % model .ne. NO_TURBULENCE_MODEL) then

    ! Variables needed for all turbulence models
    call Gpu % Vector_Update_Host(Turb % y_plus)

    ! These two belong to Field_Mod
    call Gpu % Vector_Update_Host(Flow % shear)
    call Gpu % Vector_Update_Host(Flow % vort)

  end if

  end subroutine

