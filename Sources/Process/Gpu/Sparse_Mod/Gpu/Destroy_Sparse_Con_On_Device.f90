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

  !$acc exit data delete(Acon % fc)
  !$acc exit data delete(Acon % row)
  !$acc exit data delete(Acon % col)
  !$acc exit data delete(Acon % dia)
  !$acc exit data delete(Acon % pos)

# if T_FLOWS_GPU == 1
    Gpu % gb_used = Gpu % gb_used - (  real(sizeof(Acon % fc))     &
                                     + real(sizeof(Acon % row))    &
                                     + real(sizeof(Acon % col))    &
                                     + real(sizeof(Acon % dia))    &
                                     + real(sizeof(Acon % pos))) / GIGABYTE
    print '(a,f7.3,a)', ' # '//__FILE__//' :', Gpu % gb_used, ' GB on device'
# endif

  end subroutine

