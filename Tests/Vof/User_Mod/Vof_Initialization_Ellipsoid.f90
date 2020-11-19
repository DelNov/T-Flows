!==============================================================================!
  subroutine Vof_Initialization_Ellipsoid(mult)
!------------------------------------------------------------------------------!
!   Initialize dependent variables.  (It is a bit of a mess still)             !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Work_Mod, only: prelim_vof => r_cell_01,  &
                      min_dist   => r_cell_02,  &
                      max_dist   => r_cell_03,  &
                      dist_node  => r_node_01
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Multiphase_Type), target :: mult
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: grid
  integer                   :: c, n, fu, i_nod
  integer                   :: ee, n_ellipses
  real                      :: radius_1, radius_2, radius_3
  real                      :: cent_x, cent_y, cent_z, dist_norm
!==============================================================================!

  ! First take aliases
  grid => mult % pnt_grid

  ! Open file to read Ellipsoid parameters:
  call File_Mod_Open_File_For_Reading('ellipsoid_parameters.ini', fu)

  call File_Mod_Read_Line(fu)
  read(line % tokens(1), *) n_ellipses

  do ee = 1, n_ellipses

    ! Initialize working arrays
    prelim_vof(:) = 0.0

    call File_Mod_Read_Line(fu)
    read(line % tokens(1), *) radius_1
    read(line % tokens(2), *) radius_2
    read(line % tokens(3), *) radius_3

    call File_Mod_Read_Line(fu)
    read(line % tokens(1), *) cent_x
    read(line % tokens(2), *) cent_y
    read(line % tokens(3), *) cent_z

    ! Normalized distance from ellipsoid center in nodes
    dist_node (:) = 0.0
    do n = 1, grid % n_nodes
      dist_node(n) = sqrt(  ((grid % xn(n) - cent_x) / radius_1)**2   &
                          + ((grid % yn(n) - cent_y) / radius_2)**2   &
                          + ((grid % zn(n) - cent_z) / radius_3)**2)
    end do

    ! Minimum and maximum normalized distance in cells
    min_dist(:) = +HUGE
    max_dist(:) = -HUGE
    do c = 1, grid % n_cells
      do i_nod = 1, grid % cells_n_nodes(c)
        n = grid % cells_n(i_nod, c)

        min_dist(c)= min(dist_node(n), min_dist(c))
        max_dist(c)= max(dist_node(n), max_dist(c))
      end do
    end do

    ! Simply interpolate linearly
    do c = 1, grid % n_cells

      ! Since surface is at 1.0 this checks if cell crosses the surface
      if (min_dist(c) < 1.0 .and. max_dist(c) > 1.0) then
        prelim_vof(c) = (1.0 - min_dist(c))  &
                      / (max_dist(c)-min_dist(c))

      else if (max_dist(c) <= 1.0) then
        prelim_vof(c) = 1.0
      end if

    end do

    ! Precision (what on earth?)
    do c = 1, grid % n_cells
      mult % vof % n(c) = max(prelim_vof(c), mult % vof % n(c))
    end do

  end do

  close(fu)

  end subroutine
