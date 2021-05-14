!==============================================================================!
  subroutine Vof_Exact_Cylinder(Vof,                     &
                                c,                       &
                                p1_x, p1_y, p1_z,        &
                                p2_x, p2_y, p2_z,        &
                                radius, height,          &
                                min_max_c1, min_max_c2,  &
                                vof_int)
!------------------------------------------------------------------------------!
!   Computes volume fraction of cell at interface                              !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Vof_Type), target :: Vof
  integer                :: c
  real                   :: p1_x, p1_y, p1_z
  real                   :: p2_x, p2_y, p2_z
  real                   :: radius, height
  real                   :: vof_int
  real                   :: min_max_c1, min_max_c2
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: grid
  integer                  :: nod, n_int, n_tot, n_tot_int, fu
  real                     :: xmin, xmax, ymin, ymax, zmin, zmax
  real                     :: zminprim, zmaxprim
  real                     :: area1, area2, area3, thick_y, h
  real                     :: dist
!==============================================================================!

  ! First take aliasesd
  grid => Vof % pnt_grid

  ! find bounding box:
  xmin =  HUGE; ymin =  HUGE; zmin =  HUGE;
  xmax = -HUGE; ymax = -HUGE; zmax = -HUGE;

  do nod = 1, grid % cells_n_nodes(c)
    xmin = min(xmin, grid % xn(grid % cells_n(nod,c)))
    ymin = min(ymin, grid % yn(grid % cells_n(nod,c)))
    zmin = min(zmin, grid % zn(grid % cells_n(nod,c)))
    xmax = max(xmax, grid % xn(grid % cells_n(nod,c)))
    ymax = max(ymax, grid % yn(grid % cells_n(nod,c)))
    zmax = max(zmax, grid % zn(grid % cells_n(nod,c)))
  end do

  thick_y = ymax-ymin

  ! Area with open top
  h = zmin - p1_z

  if (h >= 0.0) then   ! Cell is above cylinder center, no problem
    call Vof_Area_Square_Circle(xmin, xmax, p1_x, p1_z, radius, h, area1)

    h = zmax - p1_z

    call Vof_Area_Square_Circle(xmin, xmax, p1_x, p1_z, radius, h, area2)

    vof_int = ( area1 - area2 ) / ((zmax - zmin) * (xmax - xmin))

    if (min_max_c2 <= 1.0) then
      vof_int = 1.0
    end if

  else
    if (zmax > p1_z) then   ! Cell intersects cylinder center in z
      h = 0.0
      call Vof_Area_Square_Circle(xmin, xmax, p1_x, p1_z, radius, h, area1)

      h = zmax - p1_z

      call Vof_Area_Square_Circle(xmin, xmax, p1_x, p1_z, radius, h, area2)

      h = p1_z - zmin
      call Vof_Area_Square_Circle(xmin, xmax, p1_x, p1_z, radius, h, area3)

      vof_int = ( 2.0 * area1 - area2 - area3)/ ((zmax - zmin) * (xmax - xmin))

    else   ! Cell below cylinder center in z
      zminprim = 2.0 * p1_z - zmax
      zmaxprim = zminprim + zmax- zmin
      zmin = zminprim
      zmax = zmaxprim

      h = zmin - p1_z

      call Vof_Area_Square_Circle(xmin, xmax, p1_x, p1_z, radius, h, area1)

      h = zmax - p1_z

      call Vof_Area_Square_Circle(xmin, xmax, p1_x, p1_z, radius, h, area2)

      vof_int = ( area1 - area2 ) / ((zmax - zmin) * (xmax - xmin))

    end if
  end if

  end subroutine
