!==============================================================================!
  subroutine User_Mod_Initialize_Variables(Flow, Turb, Vof, Swarm, Sol)
!------------------------------------------------------------------------------!
!   Case-dependent initialization of VOF variable.                             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),  target :: Flow
  type(Turb_Type),   target :: Turb
  type(Vof_Type),    target :: Vof
  type(Swarm_Type),  target :: Swarm
  type(Solver_Type), target :: Sol
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: Grid
  type(Var_Type),   pointer :: fun
  real,             pointer :: dt
  integer                   :: c, n, i_nod, e, n_ellipses, fu
  real                      :: radius_x, radius_y, radius_z
  real                      :: cent_x, cent_y, cent_z, dist_norm
  real, contiguous, pointer :: prelim_vof(:), min_dist(:), max_dist(:)
  real, contiguous, pointer :: dist_node(:)
!==============================================================================!

  call Work % Connect_Real_Cell(prelim_vof, min_dist, max_dist)
  call Work % Connect_Real_Node(dist_node)

  ! Take aliases
  Grid => Flow % pnt_grid
  fun  => Vof % fun
  dt   => Flow % dt

  !---------------------------------!
  !   Initialize the VOF function   !
  !---------------------------------!

  ! Initialize the whole domain as 0.0
  fun % n(:) = 0.0

  ! Open file to read Ellipsoid parameters:
  call File % Open_For_Reading_Ascii('ellipsoid_parameters.ini', fu)

  call File % Read_Line(fu)
  read(line % tokens(1), *) n_ellipses

  do e = 1, n_ellipses

    ! Initialize working arrays
    prelim_vof(:) = 0.0

    ! Read line with radii
    call File % Read_Line(fu)
    read(line % tokens(1), *) radius_x
    read(line % tokens(2), *) radius_y
    read(line % tokens(3), *) radius_z

    ! Read line with coordinates of the elliposoid's center
    call File % Read_Line(fu)
    read(line % tokens(1), *) cent_x
    read(line % tokens(2), *) cent_y
    read(line % tokens(3), *) cent_z

    ! Normalized distance from ellipsoid center in nodes
    dist_node (:) = 0.0
    do n = 1, Grid % n_nodes
      dist_node(n) = sqrt(  ((Grid % xn(n) - cent_x) / radius_x)**2   &
                          + ((Grid % yn(n) - cent_y) / radius_y)**2   &
                          + ((Grid % zn(n) - cent_z) / radius_z)**2)
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

    ! This is useful if more elliposoids are
    ! defined and they intersect each other
    do c = 1, Grid % n_cells
      Vof % fun % n(c) = max(prelim_vof(c), Vof % fun % n(c))
    end do

  end do

  close(fu)

  ! Update buffer values
  call Grid % Exchange_Cells_Real(fun % n)

  ! Set old values to be the same as new ones
  fun % o(:) = fun % n(:)

  call Work % Disconnect_Real_Cell(prelim_vof, min_dist, max_dist)
  call Work % Disconnect_Real_Node(dist_node)

  end subroutine
