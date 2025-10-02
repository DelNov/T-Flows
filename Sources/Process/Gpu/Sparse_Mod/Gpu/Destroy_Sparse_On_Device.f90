!==============================================================================!
  subroutine Destroy_Sparse_On_Device(A)
!------------------------------------------------------------------------------!
!>  Destroys a sparse-matrix on the GPU, without copying it back to CPU.
!------------------------------------------------------------------------------!
!   Note: if you wanted to copy it before destroying, change delete to copyout !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Sparse_Type) :: A  !! parent connectivity matrix to destroy
!-----------------------[Avoid unused argument warning]------------------------!
# if T_FLOWS_GPU == 0
    Unused(A)
# endif
!==============================================================================!

  ! These entries used to belong to Sparse_Con_Type
  call Gpu % Vector_Real_Destroy_On_Device(A % fc)
  call Gpu % Vector_Int_Destroy_On_Device(A % row)
  call Gpu % Vector_Int_Destroy_On_Device(A % col)
  call Gpu % Vector_Int_Destroy_On_Device(A % dia)
  call Gpu % Matrix_Int_Destroy_On_Device(A % pos)

  ! These entries used to belong to Sparse_Val_Type
  call Gpu % Vector_Real_Destroy_On_Device(A % val)
  call Gpu % Vector_Real_Destroy_On_Device(A % d_inv)

  end subroutine

