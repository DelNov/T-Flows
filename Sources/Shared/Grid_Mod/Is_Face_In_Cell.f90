!==============================================================================!
  function Is_Face_In_Cell(Grid, s, c)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  logical             :: Is_Face_In_Cell
  class(Grid_Type)    :: Grid
  integer, intent(in) :: s, c
!-----------------------------------[Locals]-----------------------------------!
  integer :: i, m, n, cnt
!==============================================================================!

  m = Grid % faces_n_nodes(s)
  n = abs(Grid % cells_n_nodes(c))

  ! Obviously, faces_n must fit into cells_n
  Assert(m .le. n)

  ! Browse through elements of "faces_n" and every time you find a
  ! node in "faces_n"which is also in "cells_n", increase the count
  cnt = 0
  do i = 1, m
    if(any(Grid % cells_n(1:n, c) .eq. Grid % faces_n(i, s))) then
      cnt = cnt + 1
    end if
  end do

  ! If all members from "faces_n" have been found, return true
  if(cnt .eq. m) then
    Is_Face_In_Cell = .true.
  else
    Is_Face_In_Cell = .false.
  end if

  end function

