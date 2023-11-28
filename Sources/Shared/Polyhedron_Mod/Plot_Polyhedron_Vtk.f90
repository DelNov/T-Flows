!==============================================================================!
  subroutine Plot_Polyhedron_Vtk(Pol, head, rank)
!------------------------------------------------------------------------------!
!>  This subroutine is designed to visually inspect polyhedrons extracted from
!>  computational grids for further processin with Isoap library.  It takes
!>  a Polyhedron_Type object, a string head, and an integer rank as inputs.
!>  The subroutine generates a .vtk file named using head and rank, and
!>  writes the polyhedron data into this file.  This allows for a graphical
!>  inspection of extracted polyhedra in vtk-compatible visualization tools.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Polyhedron_Type) :: Pol   !! parent class
  character(*)           :: head  !! used for name creation, like a header
  integer, intent(in)    :: rank  !! used as a supplement to the file name
!-----------------------------------[Locals]-----------------------------------!
  integer           :: ip, is, iv, fu, ndata, npoly
  character(len=80) :: filename  ! don't use SL for separate compilation
!==============================================================================!

  write(filename,'(a,"-",i7.7,".vtk")') trim(head), rank

  open(newunit=fu, file=filename)
  write(fu,'(a26)')     '# vtk DataFile Version 2.0'
  write(fu,'(a6,i7.7)') 'File: ', rank
  write(fu,'(a5)')      'ASCII'
  write(fu,*)           ' '
  write(fu,'(a16)')     'DATASET POLYDATA'

  ! Write the points out
  write(fu,'(a6,i7,a6)') 'POINTS', Pol % n_nodes, ' float'
  do ip = 1, Pol % n_nodes
    write(fu,'(3f12.6)') Pol % nodes_xyz(ip,1),  &
                         Pol % nodes_xyz(ip,2),  &
                         Pol % nodes_xyz(ip,3)
  end do

  ! Count the polygons and data
  npoly = 0
  ndata = 0
  do is = 1, Pol % n_faces
    if(Pol % faces_n_nodes(is) .gt. 0) then
      npoly = npoly + 1
      ndata = ndata + Pol % faces_n_nodes(is)+1
    end if
  end do

  ! Write polygons and data out
  write(fu,'(a8,i7,i7)') 'POLYGONS', npoly, ndata
  do is = 1, Pol % n_faces
    if(Pol % faces_n_nodes(is) .gt. 0) then
      write(fu,'(i5)') Pol % faces_n_nodes(is)
      do iv = 1, Pol % faces_n_nodes(is)
        write(fu,'(i7)') Pol % faces_n(is,iv)-1
      end do
    end if
  end do

  ! Start section for point data
  write(fu,'(a10,i7)') 'POINT_DATA', Pol % n_nodes

  ! Write out the phi data
  write(fu,'(a19)')    'SCALARS phi float 1'
  write(fu,'(a20)')    'LOOKUP_TABLE default'
  do ip = 1, Pol % n_nodes
    write(fu,'(f12.6)') Pol % phi(ip)
  end do

  ! Write out node indices
  write(fu,'(a18)')    'SCALARS node int 1'
  write(fu,'(a20)')    'LOOKUP_TABLE default'
  do ip = 1, Pol % n_nodes
    write(fu,'(i7)') ip
  end do

  ! Write out the cell number (couldn't do it at cell/polygon)
  write(fu,'(a9,i7)') 'CELL_DATA', npoly

  write(fu,'(a18)')    'SCALARS cell int 1'
  write(fu,'(a20)')    'LOOKUP_TABLE default'
  do is = 1, Pol % n_faces
    if(Pol % faces_n_nodes(is) .gt. 0) then
      write(fu,'(i9)') rank
    end if
  end do

  close(fu)

  end
