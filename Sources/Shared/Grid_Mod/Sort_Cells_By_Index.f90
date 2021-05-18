!==============================================================================!
  subroutine Sort_Cells_By_Index(Grid, indx, n)
!------------------------------------------------------------------------------!
!   Sorts array of cells according to indx.                                    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type) :: Grid
  integer          :: n, indx(n)
!-----------------------------------[Locals]-----------------------------------!
  integer              :: i
  integer, allocatable :: work_1(:)
  integer, allocatable :: work_2(:,:)
!==============================================================================!

  allocate(work_1(   n))
  allocate(work_2(24,n))

  ! Sort number of nodes in a cell
  do i = 1, n
    work_1( indx(i) ) = Grid % cells_n_nodes(i)
  end do
  do i = 1, n
    Grid % cells_n_nodes(i) = work_1( i )
  end do

  ! Sort cells_n
  do i = 1, n
    work_2(1:8, indx(i)) = Grid % cells_n(1:8, i)
  end do
  do i = 1, n
    Grid % cells_n(1:8, i) = work_2(1:8, i)
  end do

  ! Sort cells_c
  do i = 1, n
    work_2(1:24, indx(i)) = Grid % cells_c(1:24, i)
  end do
  do i = 1, n
    Grid % cells_c(1:24, i) = work_2(1:24, i)
  end do

  deallocate(work_1)
  deallocate(work_2)

  end subroutine
