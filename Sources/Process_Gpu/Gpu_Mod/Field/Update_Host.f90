!==============================================================================!
  subroutine Field_Update_Host(Gpu, Flow)
!------------------------------------------------------------------------------!
!>  Copy all the field variables (velocity, pressure, temperature, ...) you
!>  might need for post-processing back to CPU
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
    O_Print '(a)', ' # Copying field to the host'
# endif

  call Gpu % Vector_Update_Host(Flow % u % n)
  call Gpu % Vector_Update_Host(Flow % v % n)
  call Gpu % Vector_Update_Host(Flow % w % n)
  call Gpu % Vector_Update_Host(Flow % p % n)
  if(Flow % heat_transfer) then
    call Gpu % Vector_Update_Host(Flow % t % n)
  end if

  end subroutine

