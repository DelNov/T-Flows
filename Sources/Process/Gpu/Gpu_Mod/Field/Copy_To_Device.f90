!==============================================================================!
  subroutine Field_Copy_To_Device(Gpu, Flow, Turb)
!------------------------------------------------------------------------------!
!>  Copy all the field variables (velocity, pressure, temperature, ...) you
!>  might need in your simulation to GPU.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Gpu_Type)          :: Gpu   !! parent class
  type(Field_Type), target :: Flow  !! field to transfer to device
  type(Turb_Type)          :: Turb  !! to check if shear and vort are needed
!-----------------------[Avoid unused argument warning]------------------------!
  integer                  :: sc
# if T_FLOWS_GPU == 0
    Unused(Gpu)
    Unused(Turb)
    Unused(Flow)
# endif
!==============================================================================!

# if T_FLOWS_GPU == 1
    O_Print '(a)', ' #-----------------------------'
    O_Print '(a)', ' # Copying field to the device'
    O_Print '(a)', ' #-----------------------------'
# endif

  ! Copy all the field variables
  call Gpu % Vector_Real_Copy_To_Device(Flow % pp % n)
  call Gpu % Vector_Real_Copy_To_Device(Flow % p % n)
  flow_pp_n => Flow % pp % n
  flow_p_n  => Flow % p  % n
  call Gpu % Vector_Real_Copy_To_Device(Flow % u % n)
  call Gpu % Vector_Real_Copy_To_Device(Flow % v % n)
  call Gpu % Vector_Real_Copy_To_Device(Flow % w % n)
  call Gpu % Vector_Real_Copy_To_Device(Flow % u % o)
  call Gpu % Vector_Real_Copy_To_Device(Flow % v % o)
  call Gpu % Vector_Real_Copy_To_Device(Flow % w % o)
  flow_u_n => Flow % u % n
  flow_u_o => Flow % u % o
  flow_v_n => Flow % v % n
  flow_v_o => Flow % v % o
  flow_w_n => Flow % w % n
  flow_w_o => Flow % w % o
  if(Flow % u % td_scheme .eq. PARABOLIC) then
    call Gpu % Vector_Real_Copy_To_Device(Flow % u % oo)
    call Gpu % Vector_Real_Copy_To_Device(Flow % v % oo)
    call Gpu % Vector_Real_Copy_To_Device(Flow % w % oo)
    flow_u_oo => Flow % u % oo
    flow_v_oo => Flow % v % oo
    flow_w_oo => Flow % w % oo
  end if
  call Gpu % Vector_Real_Copy_To_Device(Flow % v_flux % n)
  call Gpu % Vector_Real_Copy_To_Device(Flow % v_m)
  flow_v_flux_n => Flow % v_flux % n
  flow_v_m      => Flow % v_m

  if(Flow % heat_transfer) then
    call Gpu % Vector_Real_Copy_To_Device(Flow % t % n)
    call Gpu % Vector_Real_Copy_To_Device(Flow % t % o)
    flow_t_n  => Flow % t % n
    flow_t_o  => Flow % t % o
    if(Flow % t % td_scheme .eq. PARABOLIC) then
      call Gpu % Vector_Real_Copy_To_Device(Flow % t % oo)
      flow_t_oo => Flow % t % oo
    end if
  end if

  if(Flow % n_scalars .gt. 0) then
    call Gpu % Vector_Real_Copy_To_Device(Flow % scalar(1) % n)
    call Gpu % Vector_Real_Copy_To_Device(Flow % scalar(1) % o)
    flow_scalar_01_n  => Flow % scalar(1) % n
    flow_scalar_01_o  => Flow % scalar(1) % o
    if(Flow % scalar(1) % td_scheme .eq. PARABOLIC) then
      call Gpu % Vector_Real_Copy_To_Device(Flow % scalar(1) % oo)
      flow_scalar_01_oo => Flow % scalar(1) % oo
    end if
  end if


  ! You are going to need physical properties as well
  call Gpu % Vector_Real_Copy_To_Device(Flow % viscosity)
  call Gpu % Vector_Real_Copy_To_Device(Flow % density)
  flow_viscosity => Flow % viscosity
  flow_density   => Flow % density
  if(Flow % heat_transfer) then
    call Gpu % Vector_Real_Copy_To_Device(Flow % conductivity)
    call Gpu % Vector_Real_Copy_To_Device(Flow % capacity)
    flow_conductivity => Flow % conductivity
    flow_capacity     => Flow % capacity
  end if
  if(Flow % n_scalars .gt. 0) then
    call Gpu % Vector_Real_Copy_To_Device(Flow % diffusivity)
    flow_diffusivity => Flow % diffusivity
  end if

  ! Also copy whatever you need to compute gradients
  call Gpu % Matrix_Real_Copy_To_Device(Flow % grad_c2c)
  flow_grad_c2c => Flow % grad_c2c
  call Gpu % Vector_Real_Create_On_Device(Flow % phi_x)
  call Gpu % Vector_Real_Create_On_Device(Flow % phi_y)
  call Gpu % Vector_Real_Create_On_Device(Flow % phi_z)
  flow_phi_x => Flow % phi_x
  flow_phi_y => Flow % phi_y
  flow_phi_z => Flow % phi_z

  ! If you are modeling turbulence, you will need shear and vorticity
  if(Turb % model .ne. NO_TURBULENCE_MODEL) then
    call Gpu % Vector_Real_Copy_To_Device(Flow % shear)
    call Gpu % Vector_Real_Copy_To_Device(Flow % vort)
    flow_shear => Flow % shear
    flow_vort  => Flow % vort
  end if

  end subroutine

