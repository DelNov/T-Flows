!==============================================================================!
  subroutine Grid_Mod_Sort_Faces_By_Index(grid, indx, n)
!------------------------------------------------------------------------------!
!   Sorts array of cells according to indx.                                    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  integer         :: n, indx(n)
!-----------------------------------[Locals]-----------------------------------!
  integer              :: i
  integer, allocatable :: work_1(:)
  integer, allocatable :: work_2(:,:)
!==============================================================================!

  allocate(work_1(  n))
  allocate(work_2(4,n))

  ! Sort number of nodes in a cell
  do i = 1, n
    work_1( indx(i) ) = grid % faces_n_nodes(i)
  end do
  do i = 1, n
    grid % faces_n_nodes(i) = work_1( i )
  end do

  ! Sort cells_c
  do i = 1, n
    work_2(1:4, indx(i)) = grid % faces_n(1:4, i)
  end do
  do i = 1, n
    grid % faces_n(1:4, i) = work_2(1:4, i)
  end do

  end subroutine
