!==============================================================================!
  subroutine Field_Update_Host(Gpu, Flow, Turb)
!------------------------------------------------------------------------------!
!>  Copy all the field variables (velocity, pressure, temperature, ...) you
!>  might need for post-processing back to CPU
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Gpu_Type)  :: Gpu   !! parent class
  type(Turb_Type)  :: Turb  !! to check shear and vorticity
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
    O_Print '(a)', ' # Copying field to the host'
# endif

  call Gpu % Vector_Update_Host(Flow % u % n)
  call Gpu % Vector_Update_Host(Flow % v % n)
  call Gpu % Vector_Update_Host(Flow % w % n)
  call Gpu % Vector_Update_Host(Flow % p % n)
  if(Flow % heat_transfer) then
    call Gpu % Vector_Update_Host(Flow % t % n)
  end if

  ! If you are modeling turbulence, you might want shear and vorticity
  if(Turb % model .ne. NO_TURBULENCE_MODEL) then
    call Gpu % Vector_Update_Host(Flow % shear)
    call Gpu % Vector_Update_Host(Flow % vort)
  end if

  end subroutine

