!==============================================================================!
  subroutine Calculate_Cell_Inertia(Grid)
!------------------------------------------------------------------------------!
!>  Calculate the cells' inertia for a given grid.
!>  Cell inertias are calculated from from node coordinates, browsing through
!>  faces and accumulating cell inertias formed by pyramids formed from face
!>  surface to cell centers.  These inertias are needed for LES_TVM model.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type) :: Grid  !! grid under consideration
!-----------------------------------[Locals]-----------------------------------!
  integer           :: c1, c2, s, i_nod, j_nod, m
  real              :: x_cell_1, y_cell_1, z_cell_1
  real              :: x_cell_2, y_cell_2, z_cell_2
  real              :: ixx, iyy, izz, ixy, ixz, iyz
  real, allocatable :: xf(:), yf(:), zf(:)
!==============================================================================!

  ! Allocate memory for face's node coordinates
  m = size(Grid % faces_n, 1)
  allocate(xf(m))
  allocate(yf(m))
  allocate(zf(m))

  ! (Re)initialize cell inertia tensors
  Grid % ixx(:) = 0.0
  Grid % iyy(:) = 0.0
  Grid % izz(:) = 0.0
  Grid % ixy(:) = 0.0
  Grid % ixz(:) = 0.0
  Grid % iyz(:) = 0.0

  ! Perform the calculation
  do s = 1, Grid % n_faces
    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)

    do i_nod = 1, Grid % faces_n_nodes(s)  ! for all face types
      xf(i_nod) = Grid % xn(Grid % faces_n(i_nod,s))
      yf(i_nod) = Grid % yn(Grid % faces_n(i_nod,s))
      zf(i_nod) = Grid % zn(Grid % faces_n(i_nod,s))
    end do

    ! First cell
    x_cell_1 = Grid % xc(c1)
    y_cell_1 = Grid % yc(c1)
    z_cell_1 = Grid % zc(c1)

    ! Browse through faces's nodes and add tetrahedron by tetrahedron
    do i_nod = 1, Grid % faces_n_nodes(s)
      j_nod = i_nod + 1;  if(j_nod > Grid % faces_n_nodes(s)) j_nod = 1

      call Math % Tet_Inertia(Grid % xf(s), Grid % yf(s), Grid % zf(s),  &
                              xf(i_nod),    yf(i_nod),    zf(i_nod),     &
                              xf(j_nod),    yf(j_nod),    zf(j_nod),     &
                              x_cell_1,     y_cell_1,     z_cell_1,      &
                              ixx, iyy, izz, ixy, ixz, iyz,              &
                              around_node = 4)

      ! Update inertia in the first cell
      Grid % ixx(c1) = Grid % ixx(c1) + ixx
      Grid % iyy(c1) = Grid % iyy(c1) + iyy
      Grid % izz(c1) = Grid % izz(c1) + izz
      Grid % ixy(c1) = Grid % ixy(c1) + ixy
      Grid % ixz(c1) = Grid % ixz(c1) + ixz
      Grid % iyz(c1) = Grid % iyz(c1) + iyz
    end do  ! i_nod

    ! Second cell
    if(c2 > 0) then
      x_cell_2 = Grid % xc(c2) + Grid % dx(s)
      y_cell_2 = Grid % yc(c2) + Grid % dy(s)
      z_cell_2 = Grid % zc(c2) + Grid % dz(s)

      ! Browse through faces's nodes and add tetrahedron by tetrahedron to vol.
      do i_nod = 1, Grid % faces_n_nodes(s)
        j_nod = i_nod + 1;  if(j_nod > Grid % faces_n_nodes(s)) j_nod = 1

        call Math % Tet_Inertia(Grid % xf(s), Grid % yf(s), Grid % zf(s),  &
                                xf(j_nod),    yf(j_nod),    zf(j_nod),     &
                                xf(i_nod),    yf(i_nod),    zf(i_nod),     &
                                x_cell_2,     y_cell_2,     z_cell_2,      &
                                ixx, iyy, izz, ixy, ixz, iyz,              &
                                around_node = 4)

        ! Update inertia in the second cell (not complete, continue)
        Grid % ixx(c2) = Grid % ixx(c2) + ixx
        Grid % iyy(c2) = Grid % iyy(c2) + iyy
        Grid % izz(c2) = Grid % izz(c2) + izz
        Grid % ixy(c2) = Grid % ixy(c2) + ixy
        Grid % ixz(c2) = Grid % ixz(c2) + ixz
        Grid % iyz(c2) = Grid % iyz(c2) + iyz
      end do  ! i_nod
    end if

  end do
  print *, '# Cell volumes calculated !'

  end subroutine
