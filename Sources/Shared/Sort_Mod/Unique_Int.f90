!==============================================================================!
  subroutine Sort_Mod_Unique_Int(values, nu)
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
!   Sorts an array of unique values in ascending order for an input array of   !
!   unsorted integers.  The original array is overwritten, the number of       !
!   unique members is returned in argument "nu".                               !
!---------------------------------[Arguments]----------------------------------!
  integer :: values(:)
  integer :: nu
!-----------------------------------[Locals]-----------------------------------!
  integer              :: i, min_val, max_val
  integer, allocatable :: unique(:)
!------------------------------------------------------------------------------!

  nu = size(values, 1)

  allocate(unique(nu)); unique(:) = 0

  min_val = minval(values) - 1
  max_val = maxval(values)

  i = 0
  do while(min_val < max_val)
    i = i+1
    min_val = minval(values, mask = values > min_val)
    unique(i) = min_val
  end do

  nu = i
  values(1:i) = unique(1:i)

  end subroutine
