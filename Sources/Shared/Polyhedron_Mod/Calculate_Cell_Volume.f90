!==============================================================================!
  subroutine Calculate_Cell_Volume(Pol, vol)
!------------------------------------------------------------------------------!
!   Calculate cell's volume, this information is not stored in Polyhedron      !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Polyhedron_Type) :: Pol
  real, intent(out)      :: vol
!-----------------------------------[Locals]-----------------------------------!
  integer :: i_nod, j_nod, i, j, s
  real    :: dv, xc, yc, zc, xf, yf, zf
!==============================================================================!

  ! Initiaize cell's volume
  vol = 0.0

  ! Calculate cell's centroid (maybe taking it from ...
  ! ... Grid would be easier albeit less general?)
  call Pol % Calculate_Cell_Centroid(s, xc, yc, zc)

  do s = 1, Pol % n_faces

    ! Calculate face's centroid, this information is not stored in Polyhedron
    call Pol % Calculate_Face_Centroid(s, xf, yf, zf)

    ! Make a round around faces' nodes and compute volumes on the run
    dv = 0.0
    do i_nod = 1, Pol % faces_n_nodes(s)
      j_nod = i_nod + 1;  if(j_nod > Pol % faces_n_nodes(s)) j_nod = 1
      i = Pol % faces_n(s, i_nod)
      j = Pol % faces_n(s, j_nod)

      dv = dv + abs(Math % Tet_Volume(xf, yf, zf,                    &
                                      Pol % nodes_xyz(i, 1),  &
                                      Pol % nodes_xyz(i, 2),  &
                                      Pol % nodes_xyz(i, 3),  &
                                      Pol % nodes_xyz(j, 1),  &
                                      Pol % nodes_xyz(j, 2),  &
                                      Pol % nodes_xyz(j, 3),  &
                                      xc, yc, zc))
    end do
    vol = vol + dv
  end do  ! through faces

  end subroutine
