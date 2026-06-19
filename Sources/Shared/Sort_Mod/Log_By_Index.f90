!==============================================================================!
  pure subroutine Log_By_Index(Sort, n, x, indx)
!------------------------------------------------------------------------------!
!>  Sorts integer array x according to index.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Sort_Type), intent(in)    :: Sort     !! parent class
  integer,          intent(in)    :: n        !! array size
  logical,          intent(inout) :: x(n)     !! array to be sorted
  integer,          intent(in)    :: indx(n)  !! index for sorting
!-----------------------------------[Locals]-----------------------------------!
  integer              :: i        !! loop variable
  logical, allocatable :: work(:)  !! work array
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Sort)
!==============================================================================!

  allocate(work(n)); work = .false.

  do i = 1, n
    work(indx(i)) = x(i)
  end do

  do i = 1, n
    x(i) = work(i)
  end do

  deallocate(work)

  end subroutine
