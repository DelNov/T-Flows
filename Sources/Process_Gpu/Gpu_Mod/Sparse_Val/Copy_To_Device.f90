!==============================================================================!
  subroutine Sparse_Val_Copy_To_Device(Gpu, Aval)
!------------------------------------------------------------------------------!
!>  Coppies a matrix to GPU (device).  It can't copy the whole derived type,
!>  but coppies its components which are needed for accelrated calculations.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Gpu_Type)       :: Gpu   !! parent class
  type(Sparse_Val_Type) :: Aval  !! value matrix to copy
!-----------------------[Avoid unused argument warning]------------------------!
# if T_FLOWS_GPU == 0
    Unused(Gpu)
    Unused(Aval)
# endif
!==============================================================================!

  !$acc enter data copyin(Aval % val)
  !$acc enter data copyin(Aval % d_inv)
  !$acc enter data copyin(Aval % v_m)

# if T_FLOWS_GPU == 1
    Gpu % gb_used = Gpu % gb_used + (  real(sizeof(Aval % val))    &
                                     + real(sizeof(Aval % d_inv))  &
                                     + real(sizeof(Aval % v_m))) / GIGABYTE
    print '(a,f7.3,a)', ' # '//__FILE__//' :', Gpu % gb_used, ' GB on device'
# endif

  end subroutine

