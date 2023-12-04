!==============================================================================!
  subroutine Array_Real(Mem, a, i, i_inc)
!------------------------------------------------------------------------------!
!>  Enlarges a real array a to include the index i.  Optional i_inc
!>  specifies the increment to increase memory in chunks, avoiding too frequent
!>  calls to memory management procedures.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Memory_Type), intent(in)    :: Mem    !! parent class
  real, allocatable,  intent(inout) :: a(:)   !! operand array
  integer,            intent(in)    :: i      !! array index
  integer, optional,  intent(in)    :: i_inc  !! index increment
!-----------------------------------[Locals]-----------------------------------!
  real, allocatable :: temp(:)
  integer           :: new_lower1, new_upper1
  integer           :: i_increment = 0
!==============================================================================!

  ! If not allocated, allocate it with the smallest range possible
  if(.not. allocated(a)) allocate(a(i:i))

  !----------------------------------------------!
  !   If array is not large enough, enlarge it   !
  !----------------------------------------------!
  if(.not. Mem % Probe_Real_Array(a, i)) then

    ! Set up the i increment
    if(present(i_inc)) then
!     Assert(i_inc > 0)
      i_increment = i_inc
    end if

    ! Calculate new bounds
    new_lower1 = min(lbound(a, 1), i - i_increment)
    new_upper1 = max(ubound(a, 1), i + i_increment)

    ! Allocate temp array with new bounds and initialize
    allocate(temp(new_lower1:new_upper1))
    temp = 0.0

    ! Copy old data to new array
    temp(lbound(a,1):ubound(a,1)) = a

    ! Swap arrays
    call move_alloc(from=temp, to=a)

  end if

  end subroutine
