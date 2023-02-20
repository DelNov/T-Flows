!==============================================================================!
  subroutine Plot_Iso_Polygons_Vtk(Iso, head, rank)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Iso_Polygons_Type) :: Iso
  character(*)             :: head
  integer, intent(in)      :: rank
!-----------------------------------[Locals]-----------------------------------!
  integer           :: i, ip, is, iv, fu
  integer           :: npoly, ndata, ntp
  character(len=80) :: filename  ! don't use SL for separate compilation
!==============================================================================!

  write(filename,'(a,"-",i7.7,".vtk")') trim(head), rank

  open(newunit=fu, file=filename)
  write(fu,'(a26)')     '# vtk DataFile Version 2.0'
  write(fu,'(a6,i7.7)') 'File: ', rank
  write(fu,'(a5)')      'ASCII'
  write(fu,*)           ' '
  write(fu,'(a16)')     'DATASET POLYDATA'

  ! Count the vertices (points)
  ntp = 0
  do i = 1, Iso % n_polys
    ntp = ntp + Iso % polys_n_verts(i)
  end do

  ! Write the points out
  write(fu,'(a6,i7,a6)') 'POINTS', NTP, ' float'
  do ip = 1, ntp
    write(fu,'(3f12.6)') Iso % verts_xyz(ip,1),  &
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
  write(fu,'(a8,i7,i7)') 'POLYGONS', npoly, ndata
  do is = 1, Iso % n_polys
    if(Iso % polys_n_verts(is) .gt. 0) then
      write(fu,'(i7)') Iso % polys_n_verts(is)
      do iv = 1, Iso % polys_n_verts(is)
        write(fu,'(i7)') Iso % polys_v(is,iv)-1
      end do
    end if
  end do
  close(fu)

  end
