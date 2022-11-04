!==============================================================================!
  subroutine Create_Hollowedcube(Pol)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Polyhedron_Type) :: Pol
!-----------------------------------[Locals]-----------------------------------!
  integer :: is, iv, iv2
  real    :: d0, d1, d14
!==============================================================================!

  ! Unit-length cube with a half-length cubic hollow in its center:
  d0 = 0.0
  d1 = 1.0
  d14 = 0.25

  Pol % n_faces = 12
  Pol % n_nodes = 16
  Pol % faces_n_nodes(1) = 4
  Pol % faces_n(1,1) = 1
  Pol % faces_n(1,2) = 2
  Pol % faces_n(1,3) = 3
  Pol % faces_n(1,4) = 4
  Pol % faces_n_nodes(2) = 4
  Pol % faces_n(2,1) = 2
  Pol % faces_n(2,2) = 1
  Pol % faces_n(2,3) = 5
  Pol % faces_n(2,4) = 6
  Pol % faces_n_nodes(3) = 4
  Pol % faces_n(3,1) = 3
  Pol % faces_n(3,2) = 2
  Pol % faces_n(3,3) = 6
  Pol % faces_n(3,4) = 7
  Pol % faces_n_nodes(4) = 4
  Pol % faces_n(4,1) = 4
  Pol % faces_n(4,2) = 3
  Pol % faces_n(4,3) = 7
  Pol % faces_n(4,4) = 8
  Pol % faces_n_nodes(5) = 4
  Pol % faces_n(5,1) = 1
  Pol % faces_n(5,2) = 4
  Pol % faces_n(5,3) = 8
  Pol % faces_n(5,4) = 5
  Pol % faces_n_nodes(6) = 4
  Pol % faces_n(6,1) = 6
  Pol % faces_n(6,2) = 5
  Pol % faces_n(6,3) = 8
  Pol % faces_n(6,4) = 7
  do is = 1,6
     do iv = 1,4
        iv2 = 4-iv+1
        Pol % faces_n(is+6,iv) = Pol % faces_n(is,iv2)+8
     end do
     Pol % faces_n_nodes(is+6) = 4
  end do

  !       7/----------/3
  !       /|         /|
  !      / |        / |
  !    8/__|______4/  | 
  !     |  |       |  |
  !     |  /6------|--/2
  !     | /        | /
  !     |/_________|/
  !     5           1
  Pol % nodes_xyz(1,1) = d1
  Pol % nodes_xyz(1,2) = d0
  Pol % nodes_xyz(1,3) = d1
  Pol % nodes_xyz(2,1) = d1
  Pol % nodes_xyz(2,2) = d0
  Pol % nodes_xyz(2,3) = d0
  Pol % nodes_xyz(3,1) = d1
  Pol % nodes_xyz(3,2) = d1
  Pol % nodes_xyz(3,3) = d0
  Pol % nodes_xyz(4,1) = d1
  Pol % nodes_xyz(4,2) = d1
  Pol % nodes_xyz(4,3) = d1
  Pol % nodes_xyz(5,1) = d0
  Pol % nodes_xyz(5,2) = d0
  Pol % nodes_xyz(5,3) = d1
  Pol % nodes_xyz(6,1) = d0
  Pol % nodes_xyz(6,2) = d0
  Pol % nodes_xyz(6,3) = d0
  Pol % nodes_xyz(7,1) = d0
  Pol % nodes_xyz(7,2) = d1
  Pol % nodes_xyz(7,3) = d0
  Pol % nodes_xyz(8,1) = d0
  Pol % nodes_xyz(8,2) = d1
  Pol % nodes_xyz(8,3) = d1

  Pol % nodes_xyz(9,1) = d1-d14
  Pol % nodes_xyz(9,2) = d0+d14
  Pol % nodes_xyz(9,3) = d1-d14
  Pol % nodes_xyz(10,1) = d1-d14
  Pol % nodes_xyz(10,2) = d0+d14
  Pol % nodes_xyz(10,3) = d0+d14
  Pol % nodes_xyz(11,1) = d1-d14
  Pol % nodes_xyz(11,2) = d1-d14
  Pol % nodes_xyz(11,3) = d0+d14
  Pol % nodes_xyz(12,1) = d1-d14
  Pol % nodes_xyz(12,2) = d1-d14
  Pol % nodes_xyz(12,3) = d1-d14
  Pol % nodes_xyz(13,1) = d0+d14
  Pol % nodes_xyz(13,2) = d0+d14
  Pol % nodes_xyz(13,3) = d1-d14
  Pol % nodes_xyz(14,1) = d0+d14
  Pol % nodes_xyz(14,2) = d0+d14
  Pol % nodes_xyz(14,3) = d0+d14
  Pol % nodes_xyz(15,1) = d0+d14
  Pol % nodes_xyz(15,2) = d1-d14
  Pol % nodes_xyz(15,3) = d0+d14
  Pol % nodes_xyz(16,1) = d0+d14
  Pol % nodes_xyz(16,2) = d1-d14
  Pol % nodes_xyz(16,3) = d1-d14

  end
