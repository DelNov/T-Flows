!==============================================================================!
  subroutine Destroy_Sparse_Val_On_Device(Aval)
!------------------------------------------------------------------------------!
!>  Destroys a sparse-matrix on the GPU, without copying it back to CPU.
!------------------------------------------------------------------------------!
!   Note: if you wanted to copy it before destroying, change delete to copyout !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Sparse_Val_Type) :: Aval  !! value matrix to destroy
!-----------------------[Avoid unused argument warning]------------------------!
# if T_FLOWS_GPU == 0
    Unused(Aval)
# endif
!==============================================================================!

  !$acc exit data delete(Aval % val)
  !$acc exit data delete(Aval % d_inv)

# if T_FLOWS_GPU == 1
    Gpu % gb_used = Gpu % gb_used - (  real(sizeof(Aval % val))    &
                                     + real(sizeof(Aval % d_inv))) / GIGABYTE
    print '(a,f7.3,a)', ' # '//__FILE__//' :', Gpu % gb_used, ' GB on device'
# endif

  end subroutine

