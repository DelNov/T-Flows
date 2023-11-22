!==============================================================================!
  subroutine Save_Vtk_Face(Grid, s, head, rank)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type)    :: Grid
  integer, intent(in) :: s
  character(*)        :: head
  integer, intent(in) :: rank
!-----------------------------------[Locals]-----------------------------------!
  integer           :: i_nod, n, fu
  character(len=80) :: filename  ! don't use SL for separate compilation
!==============================================================================!

  write(filename,'(a,"-",i9.9,".vtk")') trim(head), rank

  open(newunit=fu, file=filename)
  write(fu,'(a26)')     '# vtk DataFile Version 2.0'
  write(fu,'(a6,i7.7)') 'File: ', s
  write(fu,'(a5)')      'ASCII'
  write(fu,*)           ' '
  write(fu,'(a16)')     'DATASET POLYDATA'

  ! Write the points out
  write(fu,'(a6,i7,a6)') 'POINTS', Grid % faces_n_nodes(s), ' float'
  do i_nod = 1, Grid % faces_n_nodes(s)
    n = Grid % faces_n(i_nod, s)
    write(fu,'(3es15.6)') Grid % xn(n), Grid % yn(n), Grid % zn(n)
  end do

  ! Write polygons and data out
  write(fu,'(a8,i7,i7)') 'POLYGONS', 1, Grid % faces_n_nodes(s) + 1
  write(fu,'(i7)') Grid % faces_n_nodes(s)
  do i_nod = 1, Grid % faces_n_nodes(s)
    write(fu,'(i7)') i_nod-1
  end do

  ! Beginning of face data
  write(fu,'(a9,i7)') 'CELL_DATA', 1

  ! Write out the face surface vectors
  write(fu,'(a25)')   'SCALARS face_surf float 3'
  write(fu,'(a20)')   'LOOKUP_TABLE default'
  write(fu,'(3es15.6)') Grid % sx(s), Grid % sy(s), Grid % sz(s)

  close(fu)

  end
