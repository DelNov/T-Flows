!==============================================================================!
  subroutine Matrix_Real(Mem, a, i, j, i_inc, j_inc)
!------------------------------------------------------------------------------!
!>  Enlarges a real matrix to include the indices  i and j.  Optional i_inc
!>  and j_inc specify the increment to increase memory in chunks, avoiding too
!>  frequent calls to memory management procedures.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Memory_Type),   intent(in)  :: Mem           !! parent class
  real, allocatable,  intent(inout) :: a(:,:)        !! operand matrix
  integer, optional,  intent(in)    :: i, j          !! matrix index
  integer, optional,  intent(in)    :: i_inc, j_inc  !! size increment
!-----------------------------------[Locals]-----------------------------------!
  real, allocatable :: temp(:,:)
  integer           :: new_lower1, new_upper1
  integer           :: new_lower2, new_upper2
  integer           :: i_local     = 1, j_local     = 1
  integer           :: i_increment = 0, j_increment = 0
!==============================================================================!

  !-------------------------------!
  !   Handle optional arguments   !
  !-------------------------------!
  if(present(i)) i_local = i
  if(present(j)) j_local = j

  ! If not allocated, allocate it with the smallest range possible
  if(.not. allocated(a)) allocate(a(i_local:i_local,j_local:j_local))

  !-----------------------------------------------!
  !   If matrix is not large enough, enlarge it   !
  !-----------------------------------------------!
  if(.not. Mem % Probe_Real_Matrix(a, i_local, j_local)) then

    ! Set up the increment in i
    if(present(i_inc)) then
!     Assert(i_inc > 0)
      i_increment = i_inc
    end if

    ! Set up the increment in j
    if(present(j_inc)) then
!     Assert(j_inc > 0)
      j_increment = j_inc
    end if

    ! Calculate new bounds
    new_lower1 = min(lbound(a, 1), i_local - i_increment)
    new_upper1 = max(ubound(a, 1), i_local + i_increment)
    new_lower2 = min(lbound(a, 2), j_local - j_increment)
    new_upper2 = max(ubound(a, 2), j_local + j_increment)

    ! Allocate temp array with new bounds and initialize
    allocate(temp(new_lower1:new_upper1, new_lower2:new_upper2))
    temp = 0.0

    ! Copy old data to new array
    temp(lbound(a,1):ubound(a,1), lbound(a,2):ubound(a,2)) = a

    ! Swap arrays
    call move_alloc(from=temp, to=a)

  end if

  end subroutine
