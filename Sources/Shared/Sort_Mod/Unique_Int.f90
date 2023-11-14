!==============================================================================!
  pure subroutine Unique_Int(Sort, values, nu)
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
!>  Sorts an array of unique values in ascending order for an input array of
!>  unsorted integers.  The original array is overwritten, the number of
!>  unique members is returned in argument "nu".
!---------------------------------[Arguments]----------------------------------!
  class(Sort_Type), intent(in)    :: Sort       !! parent class
  integer,          intent(inout) :: values(:)  !! array to be sorted
  integer,          intent(out)   :: nu         !! new size of array
!-----------------------------------[Locals]-----------------------------------!
  integer              :: i, min_val, max_val  !! helping local variables
  integer, allocatable :: unique(:)            !! uniquely sorted array
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Sort)
!==============================================================================!

  ! Set initial size of unique list ...
  nu = size(values, 1)

  ! ... can't allocate without setting it
  allocate(unique(nu)); unique(:) = 0

  min_val = minval(values(1:nu)) - 1
  max_val = maxval(values(1:nu))

  i = 0
  do while(min_val < max_val)
    i = i + 1
    min_val = minval(values(1:nu), mask = values(1:nu) > min_val)
    unique(i) = min_val
  end do

  nu = i
  values(1:i) = unique(1:i)

  end subroutine
