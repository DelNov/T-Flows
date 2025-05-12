!==============================================================================!
  subroutine Copy_Sparse_Con_To_Device(Acon)
!------------------------------------------------------------------------------!
!>  Coppies a matrix to GPU (device).  It can't copy the whole derived type,
!>  but coppies its components which are needed for accelrated calculations.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Sparse_Con_Type) :: Acon  !! parent connectivity matrix to copy
!-----------------------[Avoid unused argument warning]------------------------!
# if T_FLOWS_GPU == 0
    Unused(Acon)
# endif
!==============================================================================!

  call Gpu % Vector_Real_Copy_To_Device(Acon % fc)
  call Gpu % Vector_Int_Copy_To_Device(Acon % row)
  call Gpu % Vector_Int_Copy_To_Device(Acon % col)
  call Gpu % Vector_Int_Copy_To_Device(Acon % dia)
  call Gpu % Matrix_Int_Copy_To_Device(Acon % pos)

  end subroutine

