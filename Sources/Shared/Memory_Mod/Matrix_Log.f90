!==============================================================================!
  subroutine Matrix_Log(Mem, a, i, j, i_range, j_range, i_inc, j_inc)
!------------------------------------------------------------------------------!
!>  Enlarges a logical matrix to include the indices i and j, or ranges
!>  of indices specified in i_range and j_range.  For each dimension, either
!>  an index or a range should be specified, but it is fine to specify one
!>  dimension with index and the other with range.  Optional i_inc and j_inc
!>  specify the increment to increase memory in chunks, avoiding too frequent
!>  calls to memory management procedures.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Memory_Type),   intent(in)    :: Mem           !! parent class
  logical, allocatable, intent(inout) :: a(:,:)        !! operand matrix
  integer, optional,    intent(in)    :: i, j          !! matrix index
  integer, optional,    intent(in)    :: i_range(2)    !! matrix range in i
  integer, optional,    intent(in)    :: j_range(2)    !! matrix range in j
  integer, optional,    intent(in)    :: i_inc, j_inc  !! size increment
!-----------------------------------[Locals]-----------------------------------!
  logical, allocatable :: temp(:,:)
  integer              :: new_i_lower
  integer              :: new_i_upper
  integer              :: new_j_lower
  integer              :: new_j_upper
  integer              :: i_lower     = 1
  integer              :: i_upper     = 1
  integer              :: i_increment = 0
  integer              :: j_lower     = 1
  integer              :: j_upper     = 1
  integer              :: j_increment = 0
  integer              :: error_code       ! allocation error code
  character(DL)        :: error_message    ! allocation error message
!==============================================================================!

  !-------------------------------------------------------!
  !   Work out the ranges in i and j (rows and columns)   !
  !-------------------------------------------------------!
  call Mem % Work_Out_I_Ranges(i, i_range, i_inc, i_lower, i_upper,  &
                               __FILE__, __LINE__)
  call Mem % Work_Out_J_Ranges(j, j_range, j_inc, j_lower, j_upper,  &
                               __FILE__, __LINE__)

  !--------------------------------------------------------------------!
  !   If not allocated, allocate it with the smallest range possible   !
  !--------------------------------------------------------------------!
  if(.not. allocated(a)) then
    allocate(a(i_lower:i_upper, j_lower:j_upper),  &
             stat=error_code, errmsg=error_message)
    if(error_code .ne. 0) then
      call Message % Error(72,                                        &
         'Failed to allocate requested memory.  Message from the '//  &
         'compiler reads: "'//trim(error_message)//'". '          //  &
         'This error is critical.  Exiting!',                         &
         file=__FILE__, line=__LINE__)
    end if
    a = .false.
  end if

  !-----------------------------------------------!
  !   If matrix is not large enough, enlarge it   !
  !-----------------------------------------------!
  if(     .not. Mem % Test_Matrix_Log(a, i_lower, j_lower)  &
     .or. .not. Mem % Test_Matrix_Log(a, i_upper, j_lower)  &
     .or. .not. Mem % Test_Matrix_Log(a, i_lower, j_upper)  &
     .or. .not. Mem % Test_Matrix_Log(a, i_upper, j_upper)  ) then

    ! Set up the increment in i
    if(present(i_inc)) then
      Assert(i_inc > 0)
      i_increment = i_inc
    end if

    ! Set up the increment in j
    if(present(j_inc)) then
      Assert(j_inc > 0)
      j_increment = j_inc
    end if

    ! Calculate new bounds
    new_i_lower = min(lbound(a, 1), i_lower - i_increment)
    new_i_upper = max(ubound(a, 1), i_upper + i_increment)
    new_j_lower = min(lbound(a, 2), j_lower - j_increment)
    new_j_upper = max(ubound(a, 2), j_upper + j_increment)

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
    temp = .false.

    ! Copy old data to new array
    temp(lbound(a,1):ubound(a,1), lbound(a,2):ubound(a,2)) = a

    ! Swap arrays
    call move_alloc(from=temp, to=a)

  end if

  end subroutine
