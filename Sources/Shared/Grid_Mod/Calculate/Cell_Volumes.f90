!==============================================================================!
  subroutine Calculate_Cell_Volumes(Grid)
!------------------------------------------------------------------------------!
!>  Calculate cell volumes from node coordinates.
!------------------------------------------------------------------------------!
!   Calculate the cell volumes from node coordinates, browsing through faces   !
!   and storing partial volumes occupied in cells formed by pyramids formed    !
!   from face surface to cell centers.  These might be needed for improvements !
!   of Rhie and Chow technique in presence of body forces.                     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type) :: Grid  !! grid under consideration
!-----------------------------------[Locals]-----------------------------------!
  integer           :: c1, c2, s, i_nod, j_nod, m
  real              :: x_cell_1, y_cell_1, z_cell_1
  real              :: x_cell_2, y_cell_2, z_cell_2
  real, allocatable :: xf(:), yf(:), zf(:)
!==============================================================================!

  ! Allocate memory for face's node coordinates
  m = size(Grid % faces_n, 1)
  allocate(xf(m))
  allocate(yf(m))
  allocate(zf(m))

  ! (Re)initialize cell volumes
  Grid % vol(:) = 0.0
  Grid % dv1(:) = 0.0
  Grid % dv2(:) = 0.0

  ! Calculate volumes
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

    ! Browse through faces's nodes and add tetrahedron by tetrahedron to volume
    do i_nod = 1, Grid % faces_n_nodes(s)
      j_nod = i_nod + 1;  if(j_nod > Grid % faces_n_nodes(s)) j_nod = 1

      Grid % dv1(s) = Grid % dv1(s) +                                     &
         abs(Math % Tet_Volume(Grid % xf(s), Grid % yf(s), Grid % zf(s),  &
                               xf(i_nod),    yf(i_nod),    zf(i_nod),     &
                               xf(j_nod),    yf(j_nod),    zf(j_nod),     &
                               x_cell_1,     y_cell_1,     z_cell_1))
    end do  ! i_nod
    Grid % vol(c1) = Grid % vol(c1) + Grid % dv1(s)

    ! Second cell
    if(c2 > 0) then
      x_cell_2 = Grid % xc(c2) + Grid % dx(s)
      y_cell_2 = Grid % yc(c2) + Grid % dy(s)
      z_cell_2 = Grid % zc(c2) + Grid % dz(s)

      ! Browse through faces's nodes and add tetrahedron by tetrahedron to vol.
      do i_nod = 1, Grid % faces_n_nodes(s)
        j_nod = i_nod + 1;  if(j_nod > Grid % faces_n_nodes(s)) j_nod = 1

        Grid % dv2(s) = Grid % dv2(s) +                                     &
           abs(Math % Tet_Volume(Grid % xf(s), Grid % yf(s), Grid % zf(s),  &
                                 xf(i_nod),    yf(i_nod),    zf(i_nod),     &
                                 xf(j_nod),    yf(j_nod),    zf(j_nod),     &
                                 x_cell_2,     y_cell_2,     z_cell_2))
      end do  ! i_nod
      Grid % vol(c2) = Grid % vol(c2) + Grid % dv2(s)
    end if

  end do
  print *, '# Cell volumes calculated !'

  end subroutine
