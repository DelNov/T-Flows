!==============================================================================!
  subroutine Vector_Update_Host(Gpu, a)
!------------------------------------------------------------------------------!
!>  Copy a vector from GPU back to CPU, but do not destroy it on GPU, useful
!>  when fetching results from GPUs for saving and post-processing.
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

  call Profiler % Start('Update_Host')

  !$acc update host(a)

  call Profiler % Stop('Update_Host')

  end subroutine

