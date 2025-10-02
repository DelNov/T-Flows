!==============================================================================!
  function Is_Face_In_Cell(Grid, s, c)
!------------------------------------------------------------------------------!
!>  The function Is_Face_In_Cell is designed to determine if a face (identified
!>  by the index 's') is contained within a cell (identified by the index 'c')
!>  in a computational grid. It is primarily used for validating the integrity
!>  of cells in the grid, especially during debugging to ensure that the
!>  geometrical representation of the grid is correct.
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * Checking nodes of face against cell:                                     !
!     - The function compares the nodes making up a face (faces_n) against     !
!       the nodes of a cell (cells_n).                                         !
!   * Counting matching nodes:                                                 !
!     - It iterates through each node of the face and checks if it is also     !
!       present in the cell's nodes. For each match, a counter (cnt) is        !
!       incremented.                                                           !
!   * Determining face containment:                                            !
!     - If the number of matches (cnt) equals the number of nodes in the face  !
!       (indicating all face nodes are within the cell), it returns true;      !
!       otherwise, false.                                                      !
!   * Use in debugging:                                                        !
!     - Due to its potentially slow performance, it's recommended to use this  !
!       function primarily in debugging mode for integrity checks.             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  logical             :: Is_Face_In_Cell  !! true if the face is in the cell
  class(Grid_Type)    :: Grid             !! grid under consideration
  integer, intent(in) :: s, c             !! indices of face and cell
!-----------------------------------[Locals]-----------------------------------!
  integer :: i, m, n, cnt
!==============================================================================!

  ! Determine the number of nodes in the face and the cell
  m = Grid % faces_n_nodes(s)
  n = abs(Grid % cells_n_nodes(c))

  ! Assert the face's nodes are fewer or equal to the cell's nodes
  Assert(m .le. n)

  ! Count how many face nodes are present in the cell
  cnt = 0
  do i = 1, m
    if(any(Grid % cells_n(1:n, c) .eq. Grid % faces_n(i, s))) then
      cnt = cnt + 1
    end if
  end do

  ! Determine if all face nodes are found in the cell
  Is_Face_In_Cell = (cnt .eq. m)

  end function

