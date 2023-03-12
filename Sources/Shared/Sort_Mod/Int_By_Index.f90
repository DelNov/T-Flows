!==============================================================================!
  pure subroutine Int_By_Index(Sort, n, x, indx)
!------------------------------------------------------------------------------!
!   Sorts integer array x according to indx.                                   !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Sort_Type), intent(in)    :: Sort
  integer,          intent(in)    :: n
  integer,          intent(inout) :: x(n)
  integer,          intent(in)    :: indx(n)
!-----------------------------------[Locals]-----------------------------------!
  integer              :: i
  integer, allocatable :: work(:)
!==============================================================================!

  allocate(work(n)); work = 0

  do i = 1, n
    work(indx(i)) = x(i)
  end do

  do i = 1, n
    x(i) = work(i)
  end do

  deallocate(work)

  end subroutine
