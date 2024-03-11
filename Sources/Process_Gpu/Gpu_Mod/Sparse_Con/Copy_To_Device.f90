!==============================================================================!
  subroutine Sparse_Con_Copy_To_Device(Gpu, Acon)
!------------------------------------------------------------------------------!
!>  Coppies a matrix to GPU (device).  It can't copy the whole derived type,
!>  but coppies its components which are needed for accelrated calculations.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Gpu_Type)       :: Gpu   !! parent class
  type(Sparse_Con_Type) :: Acon  !! connectivity matrix to copy
!-----------------------[Avoid unused argument warning]------------------------!
# if T_FLOWS_GPU == 0
    Unused(Gpu)
    Unused(Acon)
# endif
!==============================================================================!

  !$acc enter data copyin(Acon % fc)
  !$acc enter data copyin(Acon % row)
  !$acc enter data copyin(Acon % col)

# if T_FLOWS_GPU == 1
    Gpu % gb_used = Gpu % gb_used + (  real(sizeof(Acon % fc))     &
                                     + real(sizeof(Acon % row))    &
                                     + real(sizeof(Acon % col))) / GIGABYTE
    print '(a,f7.3,a)', ' # '//__FILE__//' :', Gpu % gb_used, ' GB on device'
# endif

  end subroutine

