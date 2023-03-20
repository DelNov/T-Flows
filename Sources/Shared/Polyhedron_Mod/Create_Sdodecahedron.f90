!==============================================================================!
  subroutine Create_Sdodecahedron(Pol)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Polyhedron_Type) :: Pol
!------------------------------[Local parameters]------------------------------!
  integer, parameter :: NS = 200
!-----------------------------------[Locals]-----------------------------------!
  integer :: ip, ip1, ip2, ip3, is, iv, iv2, ntp0, nts0, ntsi
  real    :: a, dmod, xc, xn, xv1, xv2, yc, yn, yv1, yv2, zc, zn, zv1, zv2
  real    :: xns(NS), yns(NS), zns(NS)
!==============================================================================!

  call Pol % Create_Dodecahedron()

  do is = 1, Pol % n_faces
    ip1 = Pol % faces_n(is,1)
    ip2 = Pol % faces_n(is,2)
    ip3 = Pol % faces_n(is,3)
    xv1 = Pol % nodes_xyz(ip2,1)-Pol % nodes_xyz(ip1,1)
    yv1 = Pol % nodes_xyz(ip2,2)-Pol % nodes_xyz(ip1,2)
    zv1 = Pol % nodes_xyz(ip2,3)-Pol % nodes_xyz(ip1,3)
    xv2 = Pol % nodes_xyz(ip3,1)-Pol % nodes_xyz(ip2,1)
    yv2 = Pol % nodes_xyz(ip3,2)-Pol % nodes_xyz(ip2,2)
    zv2 = Pol % nodes_xyz(ip3,3)-Pol % nodes_xyz(ip2,3)
    xn = yv1*zv2-zv1*yv2
    yn = zv1*xv2-xv1*zv2
    zn = xv1*yv2-yv1*xv2
    dmod = (xn**2.0+yn**2.0+zn**2.0)**0.5
    xns(is) = xn/dmod
    yns(is) = yn/dmod
    zns(is) = zn/dmod
  end do

  nts0 = Pol % n_faces
  ntsi = Pol % n_faces
  ntp0 = Pol % n_nodes
  a = 0.5

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
    Pol % nodes_xyz(ntp0,1) = xc/Pol % faces_n_nodes(is)+a*XNS(is)
    Pol % nodes_xyz(ntp0,2) = yc/Pol % faces_n_nodes(is)+a*YNS(is)
    Pol % nodes_xyz(ntp0,3) = zc/Pol % faces_n_nodes(is)+a*ZNS(is)

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

  end subroutine
