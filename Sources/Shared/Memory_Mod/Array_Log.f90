!==============================================================================!
  subroutine Array_Log(Mem, a, i, i_range, i_inc)
!------------------------------------------------------------------------------!
!>  Enlarges a logical array a to include the index i, or range of indices
!>  specified by i_range.  Optional i_inc specifies the increment to increase
!>  memory in chunks, avoiding too frequent memory management procedures.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Memory_Type),   intent(in)    :: Mem         !! parent class
  logical, allocatable, intent(inout) :: a(:)        !! operand array
  integer, optional,    intent(in)    :: i           !! array index
  integer, optional,    intent(in)    :: i_range(2)  !! array range
  integer, optional,    intent(in)    :: i_inc       !! index increment
!-----------------------------------[Locals]-----------------------------------!
  logical, allocatable :: temp(:)
  integer              :: new_i_lower
  integer              :: new_i_upper
  integer              :: i_lower     = 1
  integer              :: i_upper     = 1
  integer              :: i_increment = 0
  integer              :: error_code       ! allocation error code
  character(DL)        :: error_message    ! allocation error message
!==============================================================================!

  !-------------------------------------!
  !   Work out the ranges in i (rows)   !
  !-------------------------------------!
  call Mem % Work_Out_I_Ranges(i, i_range, i_inc, i_lower, i_upper,  &
                               __FILE__, __LINE__)

  !--------------------------------------------------------------------!
  !   If not allocated, allocate it with the smallest range possible   !
  !--------------------------------------------------------------------!
  if(.not. allocated(a)) then
    allocate(a(i_lower:i_upper), stat=error_code, errmsg=error_message)
    if(error_code .ne. 0) then
      call Message % Error(72,                                        &
         'Failed to allocate requested memory.  Message from the '//  &
         'compiler reads: "'//trim(error_message)//'". '          //  &
         'This error is critical.  Exiting!',                         &
         file=__FILE__, line=__LINE__)
    end if
    a = .false.
  end if

  !----------------------------------------------!
  !   If array is not large enough, enlarge it   !
  !----------------------------------------------!
  if(     .not. Mem % Test_Array_Log(a, i_lower)  &
     .or. .not. Mem % Test_Array_Log(a, i_upper)  ) then

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
    temp = .false.

    ! Copy old data to new array
    temp(lbound(a,1):ubound(a,1)) = a

    ! Swap arrays
    call move_alloc(from=temp, to=a)

  end if

  end subroutine
