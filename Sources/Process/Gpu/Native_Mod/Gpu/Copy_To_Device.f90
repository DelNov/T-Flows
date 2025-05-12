!==============================================================================!
  subroutine Copy_Native_To_Device(Nat)
!------------------------------------------------------------------------------!
!>  Copy the preconditioner (d_inv) to the GPU and create memory for helping
!>  arrays from native solver on GPU.  As said before, It can't copy and create
!>  the whole derived type, but its components which are basic Fortran types.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Native_Type) :: Nat  !! native solver to copy to device
!-----------------------[Avoid unused argument warning]------------------------!
# if T_FLOWS_GPU == 0
    Unused(Nat)
# endif
!==============================================================================!

  ! Copy matrix and the right hand side to the device
  call Nat % C % Copy_Sparse_Con_To_Device()
  call Nat % A % Copy_Sparse_Val_To_Device()
  call Gpu % Vector_Real_Copy_To_Device(Nat % b)

  ! Create (allocate) memory for helping vectors on device
  call Gpu % Vector_Real_Create_On_Device(Nat % p)
  call Gpu % Vector_Real_Create_On_Device(Nat % q)
  call Gpu % Vector_Real_Create_On_Device(Nat % r)

  end subroutine

