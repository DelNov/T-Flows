!==============================================================================!
  subroutine Field_Copy_To_Device(Gpu, Flow)
!------------------------------------------------------------------------------!
!>  Copy all the field variables (velocity, pressure, temperature, ...) you
!>  might need in your simulation to GPU.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Gpu_Type)  :: Gpu   !! parent class
  type(Field_Type) :: Flow  !! field to transfer to device
!-----------------------[Avoid unused argument warning]------------------------!
# if T_FLOWS_GPU == 0
    Unused(Gpu)
    Unused(Flow)
# endif
!==============================================================================!

# if T_FLOWS_GPU == 0
    return
# else
    O_Print '(a)', ' #-----------------------------'
    O_Print '(a)', ' # Copying field to the device'
    O_Print '(a)', ' #-----------------------------'
# endif

  ! Copy all the field variables
  call Gpu % Vector_Real_Copy_To_Device(Flow % temp)
  call Gpu % Vector_Real_Copy_To_Device(Flow % pp % n)
  call Gpu % Vector_Real_Copy_To_Device(Flow % p % n)
  call Gpu % Vector_Real_Copy_To_Device(Flow % u % n)
  call Gpu % Vector_Real_Copy_To_Device(Flow % v % n)
  call Gpu % Vector_Real_Copy_To_Device(Flow % w % n)
  call Gpu % Vector_Real_Copy_To_Device(Flow % u % o)
  call Gpu % Vector_Real_Copy_To_Device(Flow % v % o)
  call Gpu % Vector_Real_Copy_To_Device(Flow % w % o)
  if(Flow % u % td_scheme .eq. PARABOLIC) then
    call Gpu % Vector_Real_Copy_To_Device(Flow % u % oo)
    call Gpu % Vector_Real_Copy_To_Device(Flow % v % oo)
    call Gpu % Vector_Real_Copy_To_Device(Flow % w % oo)
  end if
  call Gpu % Vector_Real_Copy_To_Device(Flow % v_flux % n)
  call Gpu % Vector_Real_Copy_To_Device(Flow % v_m)

  if(Flow % heat_transfer) then
    call Gpu % Vector_Real_Copy_To_Device(Flow % t % n)
    call Gpu % Vector_Real_Copy_To_Device(Flow % t % o)
    if(Flow % t % td_scheme .eq. PARABOLIC) then
      call Gpu % Vector_Real_Copy_To_Device(Flow % t % oo)
    end if
  end if

  ! You are going to need physical properties as well
  call Gpu % Vector_Real_Copy_To_Device(Flow % viscosity)
  call Gpu % Vector_Real_Copy_To_Device(Flow % density)
  if(Flow % heat_transfer) then
    call Gpu % Vector_Real_Copy_To_Device(Flow % conductivity)
    call Gpu % Vector_Real_Copy_To_Device(Flow % capacity)
  end if

  ! Also copy whatever you need to compute gradients
  call Gpu % Matrix_Real_Copy_To_Device(Flow % grad_c2c)
  call Gpu % Vector_Real_Create_On_Device(Flow % phi_x)
  call Gpu % Vector_Real_Create_On_Device(Flow % phi_y)
  call Gpu % Vector_Real_Create_On_Device(Flow % phi_z)

  end subroutine

