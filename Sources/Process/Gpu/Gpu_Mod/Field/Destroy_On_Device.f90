!==============================================================================!
  subroutine Field_Destroy_On_Device(Gpu, Turb, Flow)
!------------------------------------------------------------------------------!
!>  Destroys all the field variables (velocity, pressure, temperature, ...)
!>  from the device, without copying it back to the host.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Gpu_Type)  :: Gpu   !! parent class
  type(Turb_Type)  :: Turb  !! to check if shear and vort should be destroyed
  type(Field_Type) :: Flow  !! field to  to device
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
    O_Print '(a)', ' #--------------------------------'
    O_Print '(a)', ' # Destroying field on the device'
    O_Print '(a)', ' #--------------------------------'
# endif

  ! Destroy all the field variables
  call Gpu % Vector_Real_Destroy_On_Device(Flow % pp % n)
  call Gpu % Vector_Real_Destroy_On_Device(Flow % p % n)
  call Gpu % Vector_Real_Destroy_On_Device(Flow % u % n)
  call Gpu % Vector_Real_Destroy_On_Device(Flow % v % n)
  call Gpu % Vector_Real_Destroy_On_Device(Flow % w % n)
  call Gpu % Vector_Real_Destroy_On_Device(Flow % u % o)
  call Gpu % Vector_Real_Destroy_On_Device(Flow % v % o)
  call Gpu % Vector_Real_Destroy_On_Device(Flow % w % o)
  if(Flow % u % td_scheme .eq. PARABOLIC) then
    call Gpu % Vector_Real_Destroy_On_Device(Flow % u % oo)
    call Gpu % Vector_Real_Destroy_On_Device(Flow % v % oo)
    call Gpu % Vector_Real_Destroy_On_Device(Flow % w % oo)
  end if
  call Gpu % Vector_Real_Destroy_On_Device(Flow % v_flux % n)
  call Gpu % Vector_Real_Destroy_On_Device(Flow % v_m)

  if(Flow % heat_transfer) then
    call Gpu % Vector_Real_Destroy_On_Device(Flow % t % n)
    call Gpu % Vector_Real_Destroy_On_Device(Flow % t % o)
    if(Flow % t % td_scheme .eq. PARABOLIC) then
      call Gpu % Vector_Real_Destroy_On_Device(Flow % t % oo)
    end if
  end if

  ! You are not going to need physical properties any more
  call Gpu % Vector_Real_Destroy_On_Device(Flow % viscosity)
  call Gpu % Vector_Real_Destroy_On_Device(Flow % density)
  if(Flow % heat_transfer) then
    call Gpu % Vector_Real_Destroy_On_Device(Flow % conductivity)
    call Gpu % Vector_Real_Destroy_On_Device(Flow % capacity)
  end if

  ! Also destroy whatever you needed to compute gradients
  call Gpu % Matrix_Real_Destroy_On_Device(Flow % grad_c2c)
  call Gpu % Vector_Real_Destroy_On_Device(Flow % phi_x)
  call Gpu % Vector_Real_Destroy_On_Device(Flow % phi_y)
  call Gpu % Vector_Real_Destroy_On_Device(Flow % phi_z)

  ! If you are modeling turbulence, you want to destroy shear and vorticity
  if(Turb % model .ne. NO_TURBULENCE_MODEL) then
    call Gpu % Vector_Real_Destroy_On_Device(Flow % shear)
    call Gpu % Vector_Real_Destroy_On_Device(Flow % vort)
  end if

  end subroutine

