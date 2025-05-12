!==============================================================================!
  subroutine Destroy_Sparse_Con_On_Device(Acon)
!------------------------------------------------------------------------------!
!>  Destroys a sparse-matrix on the GPU, without copying it back to CPU.
!------------------------------------------------------------------------------!
!   Note: if you wanted to copy it before destroying, change delete to copyout !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Sparse_Con_Type) :: Acon  !! parent connectivity matrix to destroy
!-----------------------[Avoid unused argument warning]------------------------!
# if T_FLOWS_GPU == 0
    Unused(Acon)
# endif
!==============================================================================!

  call Gpu % Vector_Real_Destroy_On_Device(Acon % fc)
  call Gpu % Vector_Int_Destroy_On_Device(Acon % row)
  call Gpu % Vector_Int_Destroy_On_Device(Acon % col)
  call Gpu % Vector_Int_Destroy_On_Device(Acon % dia)
  call Gpu % Matrix_Int_Destroy_On_Device(Acon % pos)

  end subroutine

