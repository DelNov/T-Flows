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

  !$acc enter data copyin(Aval % val)
  !$acc enter data copyin(Aval % d_inv)

# if T_FLOWS_GPU == 1
    Gpu % gb_used = Gpu % gb_used + (  real(sizeof(Aval % val))    &
                                     + real(sizeof(Aval % d_inv))) / GIGABYTE
    print '(a,f7.3,a)', ' # '//__FILE__//' :', Gpu % gb_used, ' GB on device'
# endif

  end subroutine

