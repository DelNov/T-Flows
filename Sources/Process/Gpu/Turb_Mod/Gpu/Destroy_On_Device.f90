!==============================================================================!
  subroutine Destroy_Turb_On_Device(Turb, Flow)
!------------------------------------------------------------------------------!
!>  Destroy all the turbulence variables you don't need in GPU any more.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Turb_Type), target :: Turb  !! turbulent object to be destroyed
  type(Field_Type), target :: Flow  !! flow over which Turb is defined
!-----------------------[Avoid unused argument warning]------------------------!
# if T_FLOWS_GPU == 0
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

  !----------------------!
  !   Spalart Allmaras   !
  !----------------------!
  if(Turb % model .eq. SPALART_ALLMARAS .or.  &
     Turb % model .eq. DES_SPALART) then
    call Gpu % Vector_Real_Destroy_On_Device(Turb % vis % n)
    call Gpu % Vector_Real_Destroy_On_Device(Turb % vis % o)
    if(Turb % vis % td_scheme .eq. PARABOLIC) then
      call Gpu % Vector_Real_Destroy_On_Device(Turb % vis % oo)
    end if
    call Gpu % Vector_Real_Destroy_On_Device(Turb % vis_t)
    call Gpu % Vector_Real_Destroy_On_Device(Turb % vis_w)
    call Gpu % Vector_Real_Destroy_On_Device(Turb % h_max)
    call Gpu % Vector_Real_Destroy_On_Device(Turb % h_min)
    call Gpu % Vector_Real_Destroy_On_Device(Turb % h_w)
    if(Flow % heat_transfer) then
      call Gpu % Vector_Real_Destroy_On_Device(Turb % con_w)
    end if ! Flow % heat_transfer
  end if

  !-----------------------!
  !   Smagorinsky model   !
  !-----------------------!
  if(Turb % model .eq. LES_SMAGORINSKY) then
    call Gpu % Vector_Real_Destroy_On_Device(Turb % vis_t)
  end if

  !----------------!
  !   Wale model   !
  !----------------!
  if(Turb % model .eq. LES_WALE) then
    call Gpu % Vector_Real_Destroy_On_Device(Turb % vis_t)
    call Gpu % Vector_Real_Destroy_On_Device(Turb % wale_v)
  end if

  !------------------------------------------------!
  !   Variables needed for all turbulence models   !
  !------------------------------------------------!
  if(Turb % model .ne. NO_TURBULENCE_MODEL) then

    ! Variables needed for all turbulence models
    call Gpu % Vector_Real_Destroy_On_Device(Turb % y_plus)

    ! These two belong to Field_Mod
    call Gpu % Vector_Real_Destroy_On_Device(Flow % shear)
    call Gpu % Vector_Real_Destroy_On_Device(Flow % vort)

  end if

  end subroutine

