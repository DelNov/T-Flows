!==============================================================================!
  subroutine Vector_Real_Update_On_Host(Gpu, a, a_c_dev_ptr)
!------------------------------------------------------------------------------!
!>  Copy a vector from GPU back to CPU, but do not destroy it on GPU, useful
!>  when fetching results from GPUs for saving and post-processing.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Gpu_Type)         :: Gpu          !! parent class
  real                    :: a(:)         !! vector to copy
  type(c_ptr),   optional :: a_c_dev_ptr  !! device pointer in C-style
!-----------------------------------[Locals]-----------------------------------!
  integer :: n
!-----------------------[Avoid unused argument warning]------------------------!
# if T_FLOWS_GPU == 0
    Unused(Gpu)
    Unused(a)
# endif
!==============================================================================!

  !---------------------------------------------------------------------!
  !   If a_dev_f_ptr is not present, update data on host with OpenACC   !
  !---------------------------------------------------------------------!
  if(.not. present(a_c_dev_ptr)) then
    !$acc update host(a)

  !--------------------------------------------------------------!
  !   If a_dev_f_ptr is present, update data on host with CUDA   !
  !--------------------------------------------------------------!
  else
    n  = size(a)

    ! Copy D -> H, keeping data on device
    call cuda_copyout_double(a, a_c_dev_ptr, n)

  end if

  end subroutine

