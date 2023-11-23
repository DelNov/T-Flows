!==============================================================================!
  subroutine Sort_Faces_By_Index(Grid, indx, n)
!------------------------------------------------------------------------------!
!>  This subroutine is specialized for sorting face data in a computational
!>  grid based on a provided index array. It is primarily used in the 'Convert'
!>  phase of grid conversion to align face data according to a new ordering
!>  defined by the index array.
!------------------------------------------------------------------------------!
!   Process                                                                    !
!                                                                              !
!   * Setting Up:                                                              !
!     - Allocates temporary work arrays for sorting operations.                !
!   * Sorting Operations:                                                      !
!     - Sorts various face-related data, including the number of nodes per face,!
!       face-node connectivity (faces_n), face-cell connectivity (faces_c), and !
!       shadow faces (faces_s), according to the new index order.              !
!   * Applying Sorted Order:                                                   !
!     - Updates the face data in the grid to reflect the new sorted order.     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type) :: Grid     !! grid under consideration
  integer          :: n        !! size of the index array
  integer          :: indx(n)  !! index array
!-----------------------------------[Locals]-----------------------------------!
  integer              :: i, m
  integer, allocatable :: work_1(:)
  integer, allocatable :: work_2(:,:)
!==============================================================================!

  m = size(Grid % faces_n, 1)

  allocate(work_1   (n))
  allocate(work_2(m, n))

  ! Sort number of nodes in a face
  do i = 1, n
    work_1( indx(i) ) = Grid % faces_n_nodes(i)
  end do
  do i = 1, n
    Grid % faces_n_nodes(i) = work_1( i )
  end do

  ! Sort faces_n
  do i = 1, n
    work_2(1:m, indx(i)) = Grid % faces_n(1:m, i)
  end do
  do i = 1, n
    Grid % faces_n(1:m, i) = work_2(1:m, i)
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
