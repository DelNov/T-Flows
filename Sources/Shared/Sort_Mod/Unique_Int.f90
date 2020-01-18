!==============================================================================!
  subroutine Sort_Mod_Unique_Int(nu, values)
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
!   Sorts an array of unique values in ascending order for an input array of   !
!   unsorted integers.  The original array is overwritten, the number of       !
!   unique members is returned in argument "nu".                               !
!---------------------------------[Arguments]----------------------------------!
  integer :: nu
  integer :: values(:)
!-----------------------------------[Locals]-----------------------------------!
  integer              :: i, min_val, max_val
  integer, allocatable :: unique(:)
!------------------------------------------------------------------------------!

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
