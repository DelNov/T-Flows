!==============================================================================!
  subroutine Matrix_Int_Copy_To_Device(Gpu, a)
!------------------------------------------------------------------------------!
!>  Copy a integer matrix from CPU to GPU.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Gpu_Type) :: Gpu     !! parent class
  integer         :: a(:,:)  !! matrix to copy
!-----------------------[Avoid unused argument warning]------------------------!
# if T_FLOWS_GPU == 0
    Unused(Gpu)
    Unused(a)
# endif
!==============================================================================!

  !$acc enter data copyin(a)

# if T_FLOWS_GPU == 1
    Gpu % gb_used = Gpu % gb_used + real(sizeof(a)) / GIGABYTE
    print '(a,f7.3,a)', ' # '//__FILE__//' :', Gpu % gb_used, ' GB on device'
# endif

  end subroutine

