!==============================================================================!
  subroutine Copy_Sparse_Val_To_Device(Aval)
!------------------------------------------------------------------------------!
!>  Coppies a matrix to GPU (device).  It can't copy the whole derived type,
!>  but coppies its components which are needed for accelrated calculations.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Sparse_Val_Type) :: Aval  !! parent value matrix to copy
!-----------------------[Avoid unused argument warning]------------------------!
# if T_FLOWS_GPU == 0
    Unused(Aval)
# endif
!==============================================================================!

  call Gpu % Vector_Real_Copy_To_Device(Aval % val)
  call Gpu % Vector_Real_Copy_To_Device(Aval % d_inv)

  end subroutine

