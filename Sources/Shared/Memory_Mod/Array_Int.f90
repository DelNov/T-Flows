!==============================================================================!
  subroutine Array_Int(Mem, a, i, i_inc)
!------------------------------------------------------------------------------!
!>  Enlarges an integer array a to include the range of indices specified by i.
!>  Optional i_inc specifies the increment to increase memory in chunks,
!>  avoiding too frequent memory management procedures.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Memory_Type),   intent(in)    :: Mem    !! parent class
  integer, allocatable, intent(inout) :: a(:)   !! operand array
  integer,              intent(in)    :: i(:)   !! array range
  integer, optional,    intent(in)    :: i_inc  !! index increment
!-----------------------------------[Locals]-----------------------------------!
  integer, allocatable :: temp(:)
  integer              :: new_i_lower
  integer              :: new_i_upper
  integer              :: i_lower     = 1
  integer              :: i_upper     = 1
  integer              :: i_increment = 0
  integer              :: error_code       ! allocation error code
  character(DL)        :: error_message    ! allocation error message
!==============================================================================!

  ! Array range must be specified in an array with the length of two
  Assert(size(i) .eq. 2)

  !----------------------------------------------------------------!
  !   If not allocated, allocate, initialize and get out of here   !
  !----------------------------------------------------------------!
  if(.not. allocated(a)) then

    ! Check validity of the arguments
    Assert(i(1) .le. i(2))

    ! Allocate memory
    allocate(a(i(1):i(2)),  &
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

  !-----------------------------------------!
  !   If here, array was allocated, check   !
  !   if it is big enough and enlarge it    !
  !-----------------------------------------!
  i_lower = i(1)
  i_upper = i(2)
  if(     .not. Mem % Test_Array_Int(a, i_lower)  &
     .or. .not. Mem % Test_Array_Int(a, i_upper)  ) then

    ! Set up the increment in i
    if(present(i_inc)) then
      Assert(i_inc > 0)
      i_increment = i_inc
    end if

    ! Calculate new bounds
    new_i_lower = min(lbound(a, 1), i_lower - i_increment)
    new_i_upper = max(ubound(a, 1), i_upper + i_increment)

    ! Allocate temp array with new bounds and initialize
    allocate(temp(new_i_lower:new_i_upper),  &
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
    temp(lbound(a,1):ubound(a,1)) = a

    ! Swap arrays
    call move_alloc(from=temp, to=a)

  end if

  end subroutine
