!==============================================================================!
  subroutine Plot_Polyhedron_Vtk(Pol, head, rank)
!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Polyhedron_Type) :: Pol
  character(*)           :: head
  integer, intent(in)    :: rank
!-----------------------------------[Locals]-----------------------------------!
  integer           :: ip, is, iv, ndata, npoly
  character(len=80) :: filename  ! don't use SL for separate compilation
!==============================================================================!

  write(filename,'(a,"-",i7.7,".vtk")') trim(head), rank

  open(11, file=filename)
  write(11,'(a26)')     '# vtk DataFile Version 2.0'
  write(11,'(a6,i7.7)') 'File: ', rank
  write(11,'(a5)')      'ASCII'
  write(11,*)           ' '
  write(11,'(a16)')     'DATASET POLYDATA'

  ! Write the points out
  write(11,'(a6,i7,a6)') 'POINTS', Pol % n_nodes, ' float'
  do ip = 1, Pol % n_nodes
     write(11,'(3f12.6)') Pol % nodes_xyz(ip,1),  &
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
  write(11,'(a8,i7,i7)') 'POLYGONS', npoly, ndata
  do is = 1, Pol % n_faces
     if(Pol % faces_n_nodes(is) .gt. 0) then
        write(11,'(i7)') Pol % faces_n_nodes(is)
        do iv = 1, Pol % faces_n_nodes(is)
           write(11,'(i7)') Pol % faces_n(is,iv)-1
        end do
     end if
  end do

  ! Start section for point data
  write(11,'(a10,i7)') 'POINT_DATA', Pol % n_nodes

  ! Write out the phi data
  write(11,'(a19)')    'SCALARS phi float 1'
  write(11,'(a20)')    'LOOKUP_TABLE default'
  do ip = 1, Pol % n_nodes
     write(11,'(f12.6)') Pol % phi(ip)
  end do

  ! Write out the cell number (couldn't do it at cell/polygon)
  write(11,'(a18)')    'SCALARS cell int 1'
  write(11,'(a20)')    'LOOKUP_TABLE default'
  do ip = 1, Pol % n_nodes
     write(11,'(i9)') rank
  end do

  ! Write out node indices
  write(11,'(a18)')    'SCALARS node int 1'
  write(11,'(a20)')    'LOOKUP_TABLE default'
  do ip = 1, Pol % n_nodes
     write(11,'(i7)') ip
  end do

  close(11)

  end
