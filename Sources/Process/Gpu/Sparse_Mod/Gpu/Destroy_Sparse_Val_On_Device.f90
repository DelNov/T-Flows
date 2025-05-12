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

  call Gpu % Vector_Real_Destroy_On_Device(Aval % val)
  call Gpu % Vector_Real_Destroy_On_Device(Aval % d_inv)

  end subroutine

