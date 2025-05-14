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

  ! Copy matrix and the right hand side to the device
  call Nat % A % Destroy_Sparse_On_Device()
  call Gpu % Vector_Real_Destroy_On_Device(Nat % b)

  ! Create (allocate) memory for helping vectors on device
  call Gpu % Vector_Real_Destroy_On_Device(Nat % p)
  call Gpu % Vector_Real_Destroy_On_Device(Nat % q)
  call Gpu % Vector_Real_Destroy_On_Device(Nat % r)

  end subroutine

