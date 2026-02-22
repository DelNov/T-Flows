!==============================================================================!
  subroutine Vector_Real_Create_On_Device(Gpu, lb, ub, a)
!------------------------------------------------------------------------------!
!>  Create memory for a real vector on GPU.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Gpu_Type)     :: Gpu       !! parent class
  integer, intent(in) :: lb, ub    !! lower and upper bound
  real                :: a(lb:ub)  !! vector to create
!-----------------------[Avoid unused argument warning]------------------------!
# if T_FLOWS_GPU == 0
    Unused(Gpu)
    Unused(a)
# endif
!==============================================================================!

  !$acc enter data create(a)

# if T_FLOWS_GPU == 1
    Gpu % gb_used = Gpu % gb_used + real(sizeof(a)) / GIGABYTE
    print '(a,f7.3,a)', ' # '//__FILE__//' :', Gpu % gb_used, ' GB on device'
# endif

  end subroutine

