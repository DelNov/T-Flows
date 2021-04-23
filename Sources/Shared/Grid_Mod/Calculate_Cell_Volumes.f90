!==============================================================================!
  subroutine Grid_Mod_Calculate_Cell_Volumes(grid)
!------------------------------------------------------------------------------!
!   Calculate the cell volumes from node coordinates, browsing through faces   !
!   and storing partial volumes occupied in cells formed by pyramids formed    !
!   from face surface to cell centers.  These might be needed for improvements !
!   of Rhie and Chow technique in presence of body forces.                     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer :: c1, c2, s, i_nod, j_nod
  real    :: xf(MAX_FACES_N_NODES), yf(MAX_FACES_N_NODES), zf(MAX_FACES_N_NODES)
  real    :: x_cell_1, y_cell_1, z_cell_1
  real    :: x_cell_2, y_cell_2, z_cell_2
!==============================================================================!

  ! (Re)initialize cell volumes
  grid % vol(:) = 0.0
  grid % dv1(:) = 0.0
  grid % dv2(:) = 0.0

  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    do i_nod = 1, grid % faces_n_nodes(s)  ! for all face types
      xf(i_nod) = grid % xn(grid % faces_n(i_nod,s))
      yf(i_nod) = grid % yn(grid % faces_n(i_nod,s))
      zf(i_nod) = grid % zn(grid % faces_n(i_nod,s))
    end do

    ! First cell
    x_cell_1 = grid % xc(c1)
    y_cell_1 = grid % yc(c1)
    z_cell_1 = grid % zc(c1)

    ! Browse through faces's nodes and add tetrahedron by tetrahedron to volume
    do i_nod = 1, grid % faces_n_nodes(s)
      j_nod = i_nod + 1;  if(j_nod > grid % faces_n_nodes(s)) j_nod = 1

      grid % dv1(s) = grid % dv1(s) +                                       &
         abs(Math_Mod_Tet_Volume(grid % xf(s), grid % yf(s), grid % zf(s),  &
                                 xf(i_nod),    yf(i_nod),    zf(i_nod),     &
                                 xf(j_nod),    yf(j_nod),    zf(j_nod),     &
                                 x_cell_1,     y_cell_1,     z_cell_1))
      grid % vol(c1) = grid % vol(c1) + grid % dv1(s)
    end do  ! i_nod

    ! Second cell
    if(c2 > 0) then
      x_cell_2 = grid % xc(c2) + grid % dx(s)
      y_cell_2 = grid % yc(c2) + grid % dy(s)
      z_cell_2 = grid % zc(c2) + grid % dz(s)

      ! Browse through faces's nodes and add tetrahedron by tetrahedron to vol.
      do i_nod = 1, grid % faces_n_nodes(s)
        j_nod = i_nod + 1;  if(j_nod > grid % faces_n_nodes(s)) j_nod = 1

        grid % dv2(s) = grid % dv2(s) +                                       &
           abs(Math_Mod_Tet_Volume(grid % xf(s), grid % yf(s), grid % zf(s),  &
                                   xf(i_nod),    yf(i_nod),    zf(i_nod),     &
                                   xf(j_nod),    yf(j_nod),    zf(j_nod),     &
                                   x_cell_2,     y_cell_2,     z_cell_2))
        grid % vol(c2) = grid % vol(c2) + grid % dv2(s)
      end do  ! i_nod
    end if

  end do
  print *, '# Cell volumes calculated !'

  end subroutine
