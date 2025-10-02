!==============================================================================!
  subroutine Matrix_Int(Mem, a, i, j)
!------------------------------------------------------------------------------!
!>  Enlarges an integer matrix to include the ranges of indices specified in
!>  i and j.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Memory_Type),   intent(in)    :: Mem     !! parent class
  integer, allocatable, intent(inout) :: a(:,:)  !! operand matrix
  integer, optional,    intent(in)    :: i(:)    !! matrix range in i
  integer, optional,    intent(in)    :: j(:)    !! matrix range in j
!-----------------------------------[Locals]-----------------------------------!
  integer, allocatable :: temp(:,:)
  integer              :: new_i_lower
  integer              :: new_i_upper
  integer              :: new_j_lower
  integer              :: new_j_upper
  integer              :: i_lower     = 1
  integer              :: i_upper     = 1
  integer              :: j_lower     = 1
  integer              :: j_upper     = 1
  integer              :: error_code       ! allocation error code
  character(DL)        :: error_message    ! allocation error message
!==============================================================================!

  if(present(i)) Assert(size(i) .eq. 2)
  if(present(j)) Assert(size(j) .eq. 2)

  !----------------------------------------------------------------!
  !   If not allocated, allocate, initialize and get out of here   !
  !----------------------------------------------------------------!
  if(.not. allocated(a)) then

    ! Check validity of the arguments
    Assert(present(i) .and. present(j))
    Assert(i(1) .le. i(2))
    Assert(j(1) .le. j(2))

    ! Allocate memory
    allocate(a(i(1):i(2), j(1):j(2)),  &
             stat=error_code, errmsg=error_message)
    if(error_code .ne. 0) then
      call Message % Error(72,                                        &
         'Failed to allocate requested memory.  Message from the '//  &
         'compiler reads: "'//trim(error_message)//'". '          //  &
         'This error is critical.  Exiting!',                         &
         file=__FILE__, line=__LINE__)
    end if

    ! Initialize
    a = 0

    ! Get out
    return
  end if

  !------------------------------------------!
  !   If here, matrix was allocated, check   !
  !    if it is big enough and enlarge it    !
  !------------------------------------------!
  if(present(i)) then
    i_lower = i(1)
    i_upper = i(2)
  else
    i_lower = lbound(a,1)
    i_upper = ubound(a,1)
  end if
  if(present(j)) then
    j_lower = j(1)
    j_upper = j(2)
  else
    j_lower = lbound(a,2)
    j_upper = ubound(a,2)
  end if
  if(     .not. Mem % Test_Matrix_Int(a, i_lower, j_lower)  &
     .or. .not. Mem % Test_Matrix_Int(a, i_upper, j_lower)  &
     .or. .not. Mem % Test_Matrix_Int(a, i_lower, j_upper)  &
     .or. .not. Mem % Test_Matrix_Int(a, i_upper, j_upper)  ) then

    ! Calculate new bounds
    new_i_lower = min(lbound(a, 1), i_lower)
    new_i_upper = max(ubound(a, 1), i_upper)
    new_j_lower = min(lbound(a, 2), j_lower)
    new_j_upper = max(ubound(a, 2), j_upper)

    ! Allocate temp array with new bounds and initialize
    allocate(temp(new_i_lower:new_i_upper, new_j_lower:new_j_upper),  &
             stat=error_code, errmsg=error_message)
    if(error_code .ne. 0) then
      call Message % Error(72,                                        &
         'Failed to allocate requested memory.  Message from the '//  &
         'compiler reads: "'//trim(error_message)//'". '          //  &
         'This error is critical.  Exiting!',                         &
         file=__FILE__, line=__LINE__)
    end if
    temp = 0

    ! Copy old data to new array
    temp(lbound(a,1):ubound(a,1), lbound(a,2):ubound(a,2)) = a

    ! Swap arrays
    call move_alloc(from=temp, to=a)

  end if

  end subroutine
