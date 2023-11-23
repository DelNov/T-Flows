!==============================================================================!
  subroutine Save_Vtk_Cell(Grid, c, head, rank)
!------------------------------------------------------------------------------!
!>  This subroutine, Save_Vtk_Cell, (and her sister Save_Vtk_Face) are designed
!>  to output detailed visual representations of individual cells (or faces),
!>  from a grid, into the VTK file format. Although not presently used, these
!>  subroutines can be particularly useful in various contexts within the code
!>  in scenarios such as examining a cell where simulations starts to diverge
!>  or examining how a cell (or a face) interacts with different phases in VOF
!>  simulations with front tracking.
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * The subroutine opens a VTK file for writing and writes the headers.      !
!   * It writes the points (node coordinates) of the specified cell.           !
!   * Polygons representing the cell faces and related data (such as local     !
!     and global face numbers, surface vectors) are written.                   !
!   * The file is then closed, completing the cell visualization process.      !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type)    :: Grid  !! the grid containing the cell
  integer, intent(in) :: c     !! particular cell to be plotted
  character(*)        :: head  !! a header (title) of the file
  integer, intent(in) :: rank  !! a numerical identifier of the file (imagine
                               !! a situation in which you plot a lot of files)
!-----------------------------------[Locals]-----------------------------------!
  integer              :: fu, i_nod, l_nod, n, i_fac, s, ndata, npoly
  character(len=SL)    :: filename
  integer, allocatable :: local_node(:)
!==============================================================================!

  allocate(local_node(Grid % n_nodes))
  local_node(:) = 0

  write(filename,'(a,"-",i9.9,".vtk")') trim(head), rank

  open(newunit=fu, file=filename)
  write(fu,'(a26)')     '# vtk DataFile Version 2.0'
  write(fu,'(a6,i7.7)') 'File: ', rank
  write(fu,'(a5)')      'ASCII'
  write(fu,*)           ' '
  write(fu,'(a16)')     'DATASET POLYDATA'

  ! Write the points out
  write(fu,'(a6,i7,a6)') 'POINTS', abs(Grid % cells_n_nodes(c)), ' float'
  l_nod = 0
  do i_nod = 1, abs(Grid % cells_n_nodes(c))
    n = Grid % cells_n(i_nod, c)
    write(fu,'(3es15.6)') Grid % xn(n), Grid % yn(n), Grid % zn(n)
    l_nod = l_nod + 1
    local_node(n) = l_nod
  end do

  ! Count the polygons and data
  npoly = 0
  ndata = 0
  do i_fac = 1, Grid % cells_n_faces(c)
    s = Grid % cells_f(i_fac, c)
    npoly = npoly + 1
    ndata = ndata + Grid % faces_n_nodes(s) + 1
  end do

  ! Write polygons and data out
  write(fu,'(a8,i7,i7)') 'POLYGONS', npoly, ndata
  do i_fac = 1, Grid % cells_n_faces(c)
    s = Grid % cells_f(i_fac, c)
    write(fu,'(i5)') Grid % faces_n_nodes(s)
    do i_nod = 1, Grid % faces_n_nodes(s)
      n = Grid % faces_n(i_nod, s)
      write(fu,'(i7)') local_node(n) - 1
    end do
  end do

  ! Beginning of face data
  write(fu,'(a9,i7)') 'CELL_DATA', npoly

  ! Write out the face local number
  write(fu,'(a24)')   'SCALARS face_local int 1'
  write(fu,'(a20)')   'LOOKUP_TABLE default'
  do i_fac = 1, Grid % cells_n_faces(c)
    write(fu,'(i9)') i_fac
  end do

  ! Write out the face global number
  write(fu,'(a25)')   'SCALARS face_global int 1'
  write(fu,'(a20)')   'LOOKUP_TABLE default'
  do i_fac = 1, Grid % cells_n_faces(c)
    write(fu,'(i9)') Grid % cells_f(i_fac, c)
  end do

  ! Write out the face surface vectors
  write(fu,'(a25)')   'SCALARS face_surf float 3'
  write(fu,'(a20)')   'LOOKUP_TABLE default'
  do i_fac = 1, Grid % cells_n_faces(c)
    s = Grid % cells_f(i_fac, c)
    write(fu,'(3es15.6)') Grid % sx(s), Grid % sy(s), Grid % sz(s)
  end do

  close(fu)

  end
