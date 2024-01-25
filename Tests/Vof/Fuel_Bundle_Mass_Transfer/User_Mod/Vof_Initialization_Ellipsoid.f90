!==============================================================================!
  subroutine Vof_Initialization_Ellipsoid(Vof)
!------------------------------------------------------------------------------!
!   Initialize dependent variables.  (It is a bit of a mess still)             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Vof_Type), target :: Vof
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: Grid
  integer                   :: c, n, fu, i_nod
  integer                   :: ee, n_ellipses
  real                      :: radius_1, radius_2, radius_3
  real                      :: cent_x, cent_y, cent_z, dist_norm
  real, contiguous, pointer :: prelim_vof(:), min_dist(:), max_dist(:)
  real, contiguous, pointer :: dist_node(:)
!==============================================================================!

  call Work % Connect_Real_Cell(prelim_vof, min_dist, max_dist)
  call Work % Connect_Real_Node(dist_node)

  ! First take aliases
  Grid => Vof % pnt_grid

  ! Open file to read Ellipsoid parameters:
  call File % Open_For_Reading_Ascii('ellipsoid_parameters.ini', fu)

  call File % Read_Line(fu)
  read(line % tokens(1), *) n_ellipses

  do ee = 1, n_ellipses

    ! Initialize working arrays
    prelim_vof(:) = 0.0

    call File % Read_Line(fu)
    read(line % tokens(1), *) radius_1
    read(line % tokens(2), *) radius_2
    read(line % tokens(3), *) radius_3

    call File % Read_Line(fu)
    read(line % tokens(1), *) cent_x
    read(line % tokens(2), *) cent_y
    read(line % tokens(3), *) cent_z

    ! Normalized distance from ellipsoid center in nodes
    dist_node (:) = 0.0
    do n = 1, Grid % n_nodes
      dist_node(n) = sqrt(  ((Grid % xn(n) - cent_x) / radius_1)**2   &
                          + ((Grid % yn(n) - cent_y) / radius_2)**2   &
                          + ((Grid % zn(n) - cent_z) / radius_3)**2)
    end do

    ! Minimum and maximum normalized distance in cells
    min_dist(:) = +HUGE
    max_dist(:) = -HUGE
    do c = 1, Grid % n_cells
      do i_nod = 1, Grid % cells_n_nodes(c)
        n = Grid % cells_n(i_nod, c)

        min_dist(c)= min(dist_node(n), min_dist(c))
        max_dist(c)= max(dist_node(n), max_dist(c))
      end do
    end do

    ! Simply interpolate linearly
    do c = 1, Grid % n_cells

      ! Since surface is at 1.0 this checks if cell crosses the surface
      if (min_dist(c) < 1.0 .and. max_dist(c) > 1.0) then
        prelim_vof(c) = (1.0 - min_dist(c))  &
                      / (max_dist(c)-min_dist(c))

      else if (max_dist(c) <= 1.0) then
        prelim_vof(c) = 1.0
      end if

    end do

    ! Precision (what on earth?)
    do c = 1, Grid % n_cells
      Vof % fun % n(c) = max(prelim_vof(c), Vof % fun % n(c))
    end do

  end do

  close(fu)

  ! Reverse
  do c = 1, Grid % n_cells
    Vof % fun % n(c) = 1.0 - Vof % fun % n(c)
  end do

  call Work % Disconnect_Real_Cell(prelim_vof, min_dist, max_dist)
  call Work % Disconnect_Real_Node(dist_node)

  end subroutine
