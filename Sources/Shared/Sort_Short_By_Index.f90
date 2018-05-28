!==============================================================================!
  subroutine Sort_Short_By_Index(x, indx, n)
!------------------------------------------------------------------------------!
!   Sorts short integer array x according to indx.                             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer(kind=4) :: n, x(n), indx(n)
!-----------------------------------[Locals]-----------------------------------!
  integer(kind=4)              :: i
  integer(kind=4), allocatable :: work(:)
!==============================================================================!

  allocate(work(n)); work = 0

  do i = 1, n
    work(indx(i)) = x(i)
  end do

  do i = 1, n
    x(i)=work(i)
  end do

  deallocate(work)

  end subroutine
