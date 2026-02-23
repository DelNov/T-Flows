!==============================================================================!
  subroutine Vector_Int_Destroy_On_Device(Gpu, a, a_f_dev_ptr)
!------------------------------------------------------------------------------!
!>  Destroys an integer vector on the GPU, without copying it back to CPU.
!------------------------------------------------------------------------------!
!   Note: if you wanted to copy it before destroying, change delete to copyout !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Gpu_Type)            :: Gpu             !! parent class
  integer                    :: a(:)            !! vector to destroy
  integer, pointer, optional :: a_f_dev_ptr(:)  !! device pointer Fortran style
!-----------------------[Avoid unused argument warning]------------------------!
# if T_FLOWS_GPU == 0
    Unused(Gpu)
    Unused(a)
# endif
!==============================================================================!

  !-----------------------------------------------------------------------!
  !   If a_dev_f_ptr is not present, delete data on device with OpenACC   !
  !-----------------------------------------------------------------------!
  if(.not. present(a_f_dev_ptr)) then
    !$acc exit data delete(a)

  !----------------------------------------------------------------!
  !   If a_dev_f_ptr is present, delete data on device with CUDA   !
  !----------------------------------------------------------------!
  else
    call cuda_free(c_loc(a_f_dev_ptr))
    nullify(a_f_dev_ptr)

  end if

# if T_FLOWS_GPU == 1
    Gpu % gb_used = Gpu % gb_used - real(sizeof(a)) / GIGABYTE
    print '(a,f7.3,a)', ' # '//__FILE__//' :', Gpu % gb_used, ' GB on device'
# endif

  end subroutine

