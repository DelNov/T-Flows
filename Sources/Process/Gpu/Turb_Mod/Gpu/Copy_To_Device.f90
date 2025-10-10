!==============================================================================!
  subroutine Copy_Turb_To_Device(Turb, Flow)
!------------------------------------------------------------------------------!
!>  Copy all the turbulence variables you need in your simulation to GPU.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Turb_Type), target :: Turb  !! turbulent object to be transferrd
  type(Field_Type), target :: Flow  !! flow over which Turb is defined
!-----------------------[Avoid unused argument warning]------------------------!
# if T_FLOWS_GPU == 0
    Unused(Turb)
    Unused(Flow)
# endif
!==============================================================================!

# if T_FLOWS_GPU == 1
    O_Print '(a)', ' #-----------------------------------------'
    O_Print '(a)', ' # Copying turbulence models to the device'
    O_Print '(a)', ' #-----------------------------------------'
# endif

  !----------------------!
  !   Spalart Allmaras   !
  !----------------------!
  if(Turb % model .eq. SPALART_ALLMARAS .or.  &
     Turb % model .eq. DES_SPALART) then
    call Gpu % Vector_Real_Copy_To_Device(Turb % vis % n)
    call Gpu % Vector_Real_Copy_To_Device(Turb % vis % o)
    turb_vis_n => Turb % vis % n
    turb_vis_o => Turb % vis % o
    if(Turb % vis % td_scheme .eq. PARABOLIC) then
      call Gpu % Vector_Real_Copy_To_Device(Turb % vis % oo)
      turb_vis_oo => Turb % vis % oo
    end if
    call Gpu % Vector_Real_Copy_To_Device(Turb % vis_t)
    call Gpu % Vector_Real_Copy_To_Device(Turb % vis_w)
    call Gpu % Vector_Real_Copy_To_Device(Turb % h_max)
    call Gpu % Vector_Real_Copy_To_Device(Turb % h_min)
    call Gpu % Vector_Real_Copy_To_Device(Turb % h_w)
    turb_vis_t => Turb % vis_t
    turb_vis_w => Turb % vis_w
    turb_h_max => Turb % h_max
    turb_h_min => Turb % h_min
    turb_h_w   => Turb % h_w
  end if

  !-----------------------!
  !   Smagorinsky model   !
  !-----------------------!
  if(Turb % model .eq. LES_SMAGORINSKY) then
    call Gpu % Vector_Real_Copy_To_Device(Turb % vis_t)
    call Gpu % Vector_Real_Copy_To_Device(Turb % vis_w)
    turb_vis_t => Turb % vis_t
    turb_vis_w => Turb % vis_w
  end if

  !----------------!
  !   Wale model   !
  !----------------!
  if(Turb % model .eq. LES_WALE) then
    call Gpu % Vector_Real_Copy_To_Device(Turb % vis_t)
    call Gpu % Vector_Real_Copy_To_Device(Turb % vis_w)
    call Gpu % Vector_Real_Copy_To_Device(Turb % wale_v)
    turb_vis_t  => Turb % vis_t
    turb_vis_w  => Turb % vis_w
    turb_wale_v => Turb % wale_v
  end if

  !------------------------------------------------!
  !   Variables needed for all turbulence models   !
  !------------------------------------------------!
  if(Turb % model .ne. NO_TURBULENCE_MODEL) then

    ! Variables needed for all turbulence models
    call Gpu % Vector_Real_Copy_To_Device(Turb % y_plus)
    call Gpu % Vector_Real_Copy_To_Device(Turb % z_o)
    turb_y_plus => Turb % y_plus
    turb_z_o    => Turb % z_o

    ! These two belong to Field_Mod
    call Gpu % Vector_Real_Copy_To_Device(Flow % shear)
    call Gpu % Vector_Real_Copy_To_Device(Flow % vort)
    flow_shear => Flow % shear
    flow_vort  => Flow % vort

    ! Energy anc scalar transport
    if(Flow % heat_transfer) then
      call Gpu % Vector_Real_Copy_To_Device(Turb % con_w)
      turb_con_w => Turb % con_w
    end if ! Flow % heat_transfer
    if(Flow % n_scalars .gt. 0) then
      turb_diff_w => Turb % diff_w
    end if

  end if

  end subroutine

