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

  !$acc enter data copyin(Acon % fc)
  !$acc enter data copyin(Acon % row)
  !$acc enter data copyin(Acon % col)
  !$acc enter data copyin(Acon % dia)
  !$acc enter data copyin(Acon % pos)

# if T_FLOWS_GPU == 1
    Gpu % gb_used = Gpu % gb_used + (  real(sizeof(Acon % fc))     &
                                     + real(sizeof(Acon % row))    &
                                     + real(sizeof(Acon % col))    &
                                     + real(sizeof(Acon % dia))    &
                                     + real(sizeof(Acon % pos))) / GIGABYTE
    print '(a,f7.3,a)', ' # '//__FILE__//' :', Gpu % gb_used, ' GB on device'
# endif

  end subroutine

