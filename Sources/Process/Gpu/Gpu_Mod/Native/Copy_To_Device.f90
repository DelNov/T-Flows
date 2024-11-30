!==============================================================================!
  subroutine Native_Copy_To_Device(Gpu, Nat)
!------------------------------------------------------------------------------!
!>  Copy the preconditioner (d_inv) to the GPU and create memory for helping
!>  arrays from native solver on GPU.  As said before, It can't copy and create
!>  the whole derived type, but its components which are basic Fortran types.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Gpu_Type)   :: Gpu  !! parent class
  type(Native_Type) :: Nat  !! native solver to create
!-----------------------[Avoid unused argument warning]------------------------!
# if T_FLOWS_GPU == 0
    Unused(Gpu)
    Unused(Nat)
# endif
!==============================================================================!

  ! Copy matrix and the right hand side to the device
  call Gpu % Sparse_Con_Copy_To_Device(Nat % C)
  call Gpu % Sparse_Val_Copy_To_Device(Nat % A)
  call Gpu % Vector_Real_Copy_To_Device(Nat % b)

  ! Create (allocate) memory for helping vectors on device
  !$acc enter data create(Nat % p)
  !$acc enter data create(Nat % q)
  !$acc enter data create(Nat % r)

# if T_FLOWS_GPU == 1
    Gpu % gb_used = Gpu % gb_used + (  real(sizeof(Nat % p))      &
                                     + real(sizeof(Nat % q))      &
                                     + real(sizeof(Nat % r))) / GIGABYTE
    print '(a,f7.3,a)', ' # '//__FILE__//' :', Gpu % gb_used, ' GB on device'
# endif

  end subroutine

