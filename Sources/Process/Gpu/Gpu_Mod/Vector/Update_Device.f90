!==============================================================================!
  subroutine Vector_Update_Device(Gpu, a)
!------------------------------------------------------------------------------!
!>  Copy a vector from CPU to GPU.  It was introduced to update right hand
!>  side vectors on the device, but maybe it will find other uses.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Gpu_Type) :: Gpu   !! parent class
  real            :: a(:)  !! vector to copy
!-----------------------[Avoid unused argument warning]------------------------!
# if T_FLOWS_GPU == 0
    Unused(Gpu)
    Unused(a)
# endif
!==============================================================================!

  call Profiler % Start('Update_Device')

  !$acc update device(a)

  call Profiler % Stop('Update_Device')

  end subroutine

