!==============================================================================!
  subroutine Update_Host(Turb, Flow)
!------------------------------------------------------------------------------!
!>  Copy all the turbulence variables you need for post-processing back to CPU
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Turb_Type), target :: Turb  !! turbulent object to be update on the host
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
    O_Print '(a)', ' # Copying turbulence variables back to the host'
# endif

  !----------------------!
  !   Spalart Allmaras   !
  !----------------------!
  if(Turb % model .eq. SPALART_ALLMARAS .or.  &
     Turb % model .eq. DES_SPALART) then
    call Gpu % Vector_Update_Host(Turb % vis % n)
    call Gpu % Vector_Update_Host(Turb % vis_t)
  end if

  !-----------------------!
  !   Smagorinsky model   !
  !-----------------------!
  if(Turb % model .eq. LES_SMAGORINSKY) then
    call Gpu % Vector_Update_Host(Turb % vis_t)
  end if

  !----------------!
  !   Wale model   !
  !----------------!
  if(Turb % model .eq. LES_WALE) then
    call Gpu % Vector_Update_Host(Turb % vis_t)
    call Gpu % Vector_Update_Host(Turb % wale_v)
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

