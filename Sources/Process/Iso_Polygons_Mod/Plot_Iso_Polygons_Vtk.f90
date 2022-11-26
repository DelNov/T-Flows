!==============================================================================!
  subroutine Plot_Iso_Polygons_Vtk(Iso, ifile)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Iso_Polygons_Type) :: Iso
  integer, intent(in)      :: ifile
!-----------------------------------[Locals]-----------------------------------!
  integer           :: i, ip, is, iv
  integer           :: npoly, ndata, ntp
  character(len=80) :: filename  ! don't use SL for separate compilation
!==============================================================================!

  write(filename,'("iso-",i7.7,".vtk")') ifile
  open(11, file=filename)

  write(11,'(a26)')     '# vtk DataFile Version 2.0'
  write(11,'(a6,i7.7)') 'File: ', ifile
  write(11,'(a5)')      'ASCII'
  write(11,*)           ' '
  write(11,'(a16)')     'DATASET POLYDATA'

  ! Count the vertices (points)
  ntp = 0
  do i = 1, Iso % n_polys
     ntp = ntp + Iso % polys_n_verts(i)
  end do

  ! Write the points out
  write(11,'(a6,i7,a6)') 'POINTS', NTP, ' float'
  do ip = 1, ntp
     write(11,'(3f12.6)') Iso % verts_xyz(ip,1),  &
                          Iso % verts_xyz(ip,2),  &
                          Iso % verts_xyz(ip,3)
  end do

  ! Count the polygons and data
  npoly = 0
  ndata = 0
  do is = 1, Iso % n_polys
     if(Iso % polys_n_verts(is) .gt. 0) then
        npoly = npoly + 1
        ndata = ndata + Iso % polys_n_verts(is) + 1
     end if
  end do

  ! Write polygons and data out
  write(11,'(a8,i7,i7)') 'POLYGONS', npoly, ndata
  do is = 1, Iso % n_polys
     if(Iso % polys_n_verts(is) .gt. 0) then
        write(11,'(i7)') Iso % polys_n_verts(is)
        do iv = 1, Iso % polys_n_verts(is)
           write(11,'(i7)') Iso % polys_v(is,iv)-1
        end do
     end if
  end do
  close(11)

  end
