!==============================================================================!
  subroutine Vector_Int_Copy_To_Device(Gpu, lb, ub, a,  &
                                            a_f_dev_ptr, a_c_dev_ptr)
!------------------------------------------------------------------------------!
!>  Copy an integer vector from CPU to GPU.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Gpu_Type)            :: Gpu             !! parent class
  integer,        intent(in) :: lb, ub          !! lower and upper bound
  integer, target            :: a(lb:ub)        !! vector to copy
  integer, pointer, optional :: a_f_dev_ptr(:)  !! device pointer Fortran style
  type(c_ptr),      optional :: a_c_dev_ptr     !! device pointer in C-style
!-----------------------------------[Locals]-----------------------------------!
  integer          :: n
  logical          :: good
  integer, pointer :: a_tmp(:)     ! temporary 1:n view
!-----------------------[Avoid unused argument warning]------------------------!
# if T_FLOWS_GPU == 0
    Unused(Gpu)
    Unused(a)
# endif
!==============================================================================!

  good = (      present(a_f_dev_ptr) .and.       present(a_c_dev_ptr)) .or.  &
         (.not. present(a_f_dev_ptr) .and. .not. present(a_c_dev_ptr))
  Assert(good)

  !-----------------------------------------------------------------------!
  !   If a_dev_f_ptr is not present, create data on device with OpenACC   !
  !-----------------------------------------------------------------------!
  if(.not. present(a_f_dev_ptr)) then
    !$acc enter data copyin(a)

  !----------------------------------------------------------------!
  !   If a_dev_f_ptr is present, create data on device with CUDA   !
  !----------------------------------------------------------------!
  else
    n  = size(a)

    ! Allocate memory with CUDA helper on the device ...
    call cuda_alloc_copyin_int(a_c_dev_ptr, a, n)

    ! ... and transform the allocated pointer to Fortan style
    call c_f_pointer(a_c_dev_ptr, a_tmp, [n])

    ! Make a_f_dev_ptr point to the same memory
    ! as a_tmp, but give it indices from lb to ub
    a_f_dev_ptr(lb:ub) => a_tmp

#   if T_FLOWS_GPU == 0
      a_f_dev_ptr => a
#   endif

  end if

# if T_FLOWS_GPU == 1
    Gpu % gb_used = Gpu % gb_used + real(sizeof(a)) / GIGABYTE
    print '(a,f7.3,a)', ' # '//__FILE__//' :', Gpu % gb_used, ' GB on device'
# endif

  end subroutine

