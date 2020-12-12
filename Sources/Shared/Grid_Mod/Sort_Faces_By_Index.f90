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

  allocate(work_1                   (n))
  allocate(work_2(MAX_FACES_N_NODES, n))

  ! Sort number of nodes in a face
  do i = 1, n
    work_1( indx(i) ) = grid % faces_n_nodes(i)
  end do
  do i = 1, n
    grid % faces_n_nodes(i) = work_1( i )
  end do

  ! Sort faces_n
  do i = 1, n
    work_2(1:MAX_FACES_N_NODES, indx(i)) =  &
      grid % faces_n(1:MAX_FACES_N_NODES, i)
  end do
  do i = 1, n
    grid % faces_n(1:MAX_FACES_N_NODES, i) =  &
      work_2(1:MAX_FACES_N_NODES, i)
  end do

  ! Sort faces_c
  do i = 1, n
    work_2(1:2, indx(i)) = grid % faces_c(1:2, i)
  end do
  do i = 1, n
    grid % faces_c(1:2, i) = work_2(1:2, i)
  end do

  ! Sort faces_s
  do i = 1, n
    work_2(1, indx(i)) = grid % faces_s(i)
  end do
  do i = 1, n
    grid % faces_s(i) = work_2(1, i)
  end do

  end subroutine
