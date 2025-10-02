!==============================================================================!
  subroutine Copy_Sparse_To_Device(A)
!------------------------------------------------------------------------------!
!>  Coppies a matrix to GPU (device).  It can't copy the whole derived type,
!>  but coppies its components which are needed for accelrated calculations.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Sparse_Type) :: A  !! parent matrix to copy
!-----------------------[Avoid unused argument warning]------------------------!
# if T_FLOWS_GPU == 0
    Unused(A)
# endif
!==============================================================================!

  ! These entries used to belong to Sparse_Con_Type
  call Gpu % Vector_Real_Copy_To_Device(A % fc)
  call Gpu % Vector_Int_Copy_To_Device(A % row)
  call Gpu % Vector_Int_Copy_To_Device(A % col)
  call Gpu % Vector_Int_Copy_To_Device(A % dia)
  call Gpu % Matrix_Int_Copy_To_Device(A % pos)

    ! These entries used to belong to Sparse_Val_Type
  call Gpu % Vector_Real_Copy_To_Device(A % val)
  call Gpu % Vector_Real_Copy_To_Device(A % d_inv)

  end subroutine

