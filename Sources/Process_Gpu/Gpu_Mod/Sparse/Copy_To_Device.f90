!==============================================================================!
  subroutine Sparse_Copy_To_Device(Gpu, A)
!------------------------------------------------------------------------------!
!>  Coppies a matrix to GPU (device).  It can't copy the whole derived type,
!>  but coppies its components which are needed for accelrated calculations.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Gpu_Type)   :: Gpu  !! parent class
  type(Sparse_Type) :: A    !! matrix to copy
!-----------------------[Avoid unused argument warning]------------------------!
# if T_FLOWS_GPU == 0
    Unused(Gpu)
    Unused(A)
# endif
!==============================================================================!

  !$acc enter data copyin(A % val)
  !$acc enter data copyin(A % fc)
  !$acc enter data copyin(A % row)
  !$acc enter data copyin(A % col)
  !$acc enter data copyin(A % d_inv)
  !$acc enter data copyin(A % v_m)

# if T_FLOWS_GPU == 1
    Gpu % gb_used = Gpu % gb_used + (  real(sizeof(A % val))    &
                                     + real(sizeof(A % fc))     &
                                     + real(sizeof(A % row))    &
                                     + real(sizeof(A % col))    &
                                     + real(sizeof(A % d_inv))  &
                                     + real(sizeof(A % v_m))) / GIGABYTE
    print '(a,f7.3,a)', ' # '//__FILE__//' :', Gpu % gb_used, ' GB on device'
# endif

  end subroutine

