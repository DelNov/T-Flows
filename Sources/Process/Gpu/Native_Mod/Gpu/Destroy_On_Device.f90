!==============================================================================!
  subroutine Destroy_Native_On_Device(Nat)
!------------------------------------------------------------------------------!
!>  Destroys a native solver on the GPU, without copying it back to CPU.
!------------------------------------------------------------------------------!
!   Note: if you wanted to copy it before destroying, change delete to copyout !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Native_Type) :: Nat  !! native solver to destroy
!-----------------------[Avoid unused argument warning]------------------------!
# if T_FLOWS_GPU == 0
    Unused(Nat)
# endif
!==============================================================================!

  !$acc exit data delete(Nat % r)
  !$acc exit data delete(Nat % p)
  !$acc exit data delete(Nat % q)

# if T_FLOWS_GPU == 1
    Gpu % gb_used = Gpu % gb_used - (  real(sizeof(Nat % p))      &
                                     + real(sizeof(Nat % q))      &
                                     + real(sizeof(Nat % r))) / GIGABYTE
    print '(a,f7.3,a)', ' # '//__FILE__//' :', Gpu % gb_used, ' GB on device'
# endif

  end subroutine

