!==============================================================================!
  subroutine Create_Zigzagcell(Pol)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Polyhedron_Type) :: Pol
!-----------------------------------[Locals]-----------------------------------!
  integer :: iv, nzigs
  real    :: d0, d1, doff
!==============================================================================!

  ! Change nzigs to modify the number of zig-zag sections
  nzigs = 5
  doff = 0.1
  d0 = 0.0
  d1 = 1.0
  Pol % n_nodes = 4*nzigs
  Pol % n_faces = 2*nzigs+2
  do iv = 1,nzigs
    Pol % nodes_xyz(iv,1) = d1*(iv-1)
    Pol % nodes_xyz(iv,2) = doff+d1*mod(iv-1,2)
    Pol % nodes_xyz(iv,3) = d0
    Pol % nodes_xyz(iv+nzigs,1) = d1*(nzigs-(iv-1)-1)
    Pol % nodes_xyz(iv+nzigs,2) = d1*mod(nzigs-(iv-1)-1,2)
    Pol % nodes_xyz(iv+nzigs,3) = d0
    Pol % nodes_xyz(iv+2*nzigs,1) = d1*(iv-1)
    Pol % nodes_xyz(iv+2*nzigs,2) = doff+d1*mod(iv-1,2)
    Pol % nodes_xyz(iv+2*nzigs,3) = d1
    Pol % nodes_xyz(iv+3*nzigs,1) = d1*(nzigs-(iv-1)-1)
    Pol % nodes_xyz(iv+3*nzigs,2) = d1*mod(nzigs-(iv-1)-1,2)
    Pol % nodes_xyz(iv+3*nzigs,3) = d1
  end do
  Pol % faces_n_nodes(1) = 2*nzigs
  Pol % faces_n_nodes(2) = 2*nzigs
  do iv = 1,2*nzigs
    Pol % faces_n(1,iv) = iv
    Pol % faces_n(2,iv) = 4*nzigs-(iv-1)
  end do
  iv = 1
  Pol % faces_n_nodes(3) = 4
  Pol % faces_n(3,1) = iv
  Pol % faces_n(3,2) = 2*nzigs-(iv-1)
  Pol % faces_n(3,3) = 4*nzigs-(iv-1)
  Pol % faces_n(3,4) = iv+2*nzigs
  iv = nzigs
  Pol % faces_n_nodes(4) = 4
  Pol % faces_n(4,4) = iv
  Pol % faces_n(4,3) = 2*nzigs-(iv-1)
  Pol % faces_n(4,2) = 4*nzigs-(iv-1)
  Pol % faces_n(4,1) = iv+2*nzigs
  do iv = 1,nzigs-1
    Pol % faces_n_nodes(iv+4) = 4
    Pol % faces_n(iv+4,1) = iv
    Pol % faces_n(iv+4,2) = iv+2*nzigs
    Pol % faces_n(iv+4,3) = iv+2*nzigs+1
    Pol % faces_n(iv+4,4) = iv+1
    Pol % faces_n_nodes(iv+4+(nzigs-1)) = 4
    Pol % faces_n(iv+4+(nzigs-1),1) = 2*nzigs-(iv-1)
    Pol % faces_n(iv+4+(nzigs-1),2) = 2*nzigs-(iv-1)-1
    Pol % faces_n(iv+4+(nzigs-1),3) = 2*nzigs-(iv-1)-1+2*nzigs
    Pol % faces_n(iv+4+(nzigs-1),4) = 2*nzigs-(iv-1)+2*nzigs
  end do

  end subroutine
