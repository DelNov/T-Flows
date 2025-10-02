!==============================================================================!
  subroutine Update_Field_On_Host(Flow)
!------------------------------------------------------------------------------!
!>  Copy all the field variables (velocity, pressure, temperature, ...) you
!>  might need for post-processing back to CPU
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Field_Type) :: Flow  !! parent field to transfer to host
!--------------------------------[Locals]--------------------------------------!
  integer :: sc
!-----------------------[Avoid unused argument warning]------------------------!
# if T_FLOWS_GPU == 0
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

  do sc = 1, Flow % n_scalars
    call Gpu % Vector_Update_Host(Flow % scalar(sc) % n)
  end do

  end subroutine

