!==============================================================================!
  subroutine Create_Scube(Pol)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Polyhedron_Type) :: Pol
!------------------------------[Local parameters]------------------------------!
  integer, parameter :: NS = 200
!-----------------------------------[Locals]-----------------------------------!
  integer :: ip, is, iv, iv2, ntp0, nts0, ntsi
  real    :: a, d0, d1, xc, yc, zc
  real    :: xns(NS), yns(NS), zns(NS)
!==============================================================================!

  call Pol % Create_Cube()
  d0 = 0.0
  d1 = 1.0

  xns(1) = d1
  yns(1) = d0
  zns(1) = d0
  xns(2) = d0
  yns(2) = -d1
  zns(2) = d0
  xns(3) = d0
  yns(3) = d0
  zns(3) = -d1
  xns(4) = d0
  yns(4) = d1
  zns(4) = d0
  xns(5) = d0
  yns(5) = d0
  zns(5) = d1
  xns(6) = -d1
  yns(6) = d0
  zns(6) = d0

  nts0 = Pol % n_faces
  ntsi = Pol % n_faces
  ntp0 = Pol % n_nodes
  a = 1.0

  ! Face centroid
  do is = 1,Pol % n_faces
     xc = 0.0
     yc = 0.0
     zc = 0.0
     do iv = 1,Pol % faces_n_nodes(is)
        ip = Pol % faces_n(is,iv)
        xc = xc+Pol % nodes_xyz(ip,1)
        yc = yc+Pol % nodes_xyz(ip,2)
        zc = zc+Pol % nodes_xyz(ip,3)
     end do
     ntp0 = ntp0+1
     Pol % nodes_xyz(ntp0,1) = xc/Pol % faces_n_nodes(is)+a*xns(is)
     Pol % nodes_xyz(ntp0,2) = yc/Pol % faces_n_nodes(is)+a*yns(is)
     Pol % nodes_xyz(ntp0,3) = zc/Pol % faces_n_nodes(is)+a*zns(is)

     do iv = 1,Pol % faces_n_nodes(is)
        iv2 = iv+1
        if(iv .eq. Pol % faces_n_nodes(is)) iv2 = 1
        nts0 = nts0+1
        Pol % faces_n_nodes(nts0) = 3
        Pol % faces_n(nts0,1) = Pol % faces_n(is,iv)
        Pol % faces_n(nts0,2) = Pol % faces_n(is,iv2)
        Pol % faces_n(nts0,3) = ntp0
     end do
  end do

  Pol % n_faces = nts0-ntsi
  do is = 1,Pol % n_faces
     Pol % faces_n_nodes(is) = Pol % faces_n_nodes(ntsi+is)
     do iv = 1,Pol % faces_n_nodes(is)
        Pol % faces_n(is,iv) = Pol % faces_n(ntsi+is,iv)
     end do
  end do
  Pol % n_nodes = ntp0

  end
