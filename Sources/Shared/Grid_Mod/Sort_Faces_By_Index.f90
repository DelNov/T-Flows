!==============================================================================!
  subroutine Sort_Faces_By_Index(Grid, indx, n)
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

  allocate(work_1                   (n))
  allocate(work_2(MAX_FACES_N_NODES, n))

  ! Sort number of nodes in a face
  do i = 1, n
    work_1( indx(i) ) = Grid % faces_n_nodes(i)
  end do
  do i = 1, n
    Grid % faces_n_nodes(i) = work_1( i )
  end do

  ! Sort faces_n
  do i = 1, n
    work_2(1:MAX_FACES_N_NODES, indx(i)) =  &
      Grid % faces_n(1:MAX_FACES_N_NODES, i)
  end do
  do i = 1, n
    Grid % faces_n(1:MAX_FACES_N_NODES, i) =  &
      work_2(1:MAX_FACES_N_NODES, i)
  end do

  ! Sort faces_c
  do i = 1, n
    work_2(1:2, indx(i)) = Grid % faces_c(1:2, i)
  end do
  do i = 1, n
    Grid % faces_c(1:2, i) = work_2(1:2, i)
  end do

  ! Sort faces_s
  do i = 1, n
    work_2(1, indx(i)) = Grid % faces_s(i)
  end do
  do i = 1, n
    Grid % faces_s(i) = work_2(1, i)
  end do

  end subroutine
