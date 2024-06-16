!==============================================================================!
  subroutine Turb_Destroy_On_Device(Gpu, Turb, Flow)
!------------------------------------------------------------------------------!
!>  Destroy all the turbulence variables you don't need in GPU any more.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Gpu_Type)  :: Gpu   !! parent class
  type(Turb_Type)  :: Turb  !! to check if wall distance should be destroyed
  type(Field_Type) :: Flow  !! grid to destroy on device
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
    O_Print '(a)', ' #-----------------------------------------------'
    O_Print '(a)', ' # Destroying turbulence variables on the device'
    O_Print '(a)', ' #-----------------------------------------------'
# endif

  !-----------------------!
  !   Smagorinsky model   !
  !-----------------------!
  if(Turb % model .eq. LES_SMAGORINSKY) then
    call Gpu % Vector_Real_Destroy_On_Device(Turb % vis_t)
  end if

  !------------------------------------------------!
  !   Variables needed for all turbulence models   !
  !------------------------------------------------!
  if(Turb % model .ne. NO_TURBULENCE_MODEL) then

    call Gpu % Vector_Real_Destroy_On_Device(Turb % y_plus)

    ! These two belong to Field_Mod
    call Gpu % Vector_Real_Destroy_On_Device(Flow % shear)
    call Gpu % Vector_Real_Destroy_On_Device(Flow % vort)

  end if

  end subroutine

