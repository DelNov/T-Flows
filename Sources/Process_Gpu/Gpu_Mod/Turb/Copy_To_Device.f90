!==============================================================================!
  subroutine Turb_Copy_To_Device(Gpu, Turb, Flow)
!------------------------------------------------------------------------------!
!>  Copy all the turbulence variables you need in your simulation to GPU.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Gpu_Type)          :: Gpu   !! parent class
  type(Turb_Type),  target :: Turb  !! turbulent object to be transferrd
  type(Field_Type), target :: Flow  !! grid to transfer to device
!-----------------------[Avoid unused argument warning]------------------------!
# if T_FLOWS_GPU == 0
    Unused(Gpu)
    Unused(Turb)
    Unused(Flow)
# endif
!==============================================================================!

# if T_FLOWS_GPU == 1
    O_Print '(a)', ' #-----------------------------------------'
    O_Print '(a)', ' # Copying turbulence models to the device'
    O_Print '(a)', ' #-----------------------------------------'
# endif

  !-----------------------!
  !   Smagorinsky model   !
  !-----------------------!
  if(Turb % model .eq. LES_SMAGORINSKY) then
    call Gpu % Vector_Real_Copy_To_Device(Turb % vis_t)
    turb_vis_t => Turb % vis_t
  end if

  !----------------!
  !   Wale model   !
  !----------------!
  if(Turb % model .eq. LES_WALE) then
    call Gpu % Vector_Real_Copy_To_Device(Turb % vis_t)
    call Gpu % Vector_Real_Copy_To_Device(Turb % wale_v)
    turb_vis_t  => Turb % vis_t
    turb_wale_v => Turb % wale_v
  end if

  !------------------------------------------------!
  !   Variables needed for all turbulence models   !
  !------------------------------------------------!
  if(Turb % model .ne. NO_TURBULENCE_MODEL) then

    ! Variables needed for all turbulence models
    call Gpu % Vector_Real_Copy_To_Device(Turb % y_plus)
    turb_y_plus => Turb % y_plus

    ! These two belong to Field_Mod
    call Gpu % Vector_Real_Copy_To_Device(Flow % shear)
    call Gpu % Vector_Real_Copy_To_Device(Flow % vort)
    flow_shear => Flow % shear
    flow_vort  => Flow % vort

  end if

  end subroutine

