!==============================================================================!
  subroutine Calculate_Cell_Volume(Polyhedron, vol)
!------------------------------------------------------------------------------!
!   Calculate cell's volume, this information is not stored in Polyhedron      !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Polyhedron_Type) :: Polyhedron
  real, intent(out)      :: vol
!-----------------------------------[Locals]-----------------------------------!
  integer :: i_nod, j_nod, i, j, s
  real    :: dv, xc, yc, zc, xf, yf, zf
!==============================================================================!

  ! Initiaize cell's volume
  vol = 0.0

  ! Calculate cell's centroid (maybe taking it from ...
  ! ... Grid would be easier albeit less general?)
  call Polyhedron % Calculate_Cell_Centroid(s, xc, yc, zc)

  do s = 1, Polyhedron % n_faces

    ! Calculate face's centroid, this information is not stored in Polyhedron
    call Polyhedron % Calculate_Face_Centroid(s, xf, yf, zf)

    ! Make a round around faces' nodes and compute volumes on the run
    dv = 0.0
    do i_nod = 1, Polyhedron % faces_n_nodes(s)
      j_nod = i_nod + 1;  if(j_nod > Polyhedron % faces_n_nodes(s)) j_nod = 1
      i = Polyhedron % faces_n(s, i_nod)
      j = Polyhedron % faces_n(s, j_nod)

      dv = dv + abs(Math % Tet_Volume(xf, yf, zf,                    &
                                      Polyhedron % nodes_xyz(i, 1),  &
                                      Polyhedron % nodes_xyz(i, 2),  &
                                      Polyhedron % nodes_xyz(i, 3),  &
                                      Polyhedron % nodes_xyz(j, 1),  &
                                      Polyhedron % nodes_xyz(j, 2),  &
                                      Polyhedron % nodes_xyz(j, 3),  &
                                      xc, yc, zc))
    end do
    vol = vol + dv
  end do  ! through faces

  end subroutine
