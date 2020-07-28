!==============================================================================!
  subroutine Vof_Interface_Ellipsoid(mult,                          &
                                    c,                              &
                                    cent_x, cent_y, cent_z,         &
                                    radius_1, radius_2, radius_3,   &
                                    vof_int)
!------------------------------------------------------------------------------!
!                Computes volume fraction of cell at interface                 !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Multiphase_Type), target :: mult
  integer                       :: c
  real                          :: cent_x, cent_y, cent_z
  real                          :: radius_1, radius_2, radius_3
  real                          :: vof_int
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),      pointer :: grid
  logical                       :: i_cell
  integer, parameter            :: n = 10000
  integer                       :: nod, n_int, n_tot, fu
  integer                       :: ee, n_cylinders, i_vari, n_vari
  real                          :: r_num, res_func
  real                          :: xmin, xmax, ymin, ymax, zmin, zmax
  real                          :: p(3)
  real                          :: vof_0, vof_tol1, vof_tol2
  real                          :: avg_x, avg_y, avg_z
  real                          :: var_comb, var_comb_0, dist_cent
  real                          :: mean_x, mean_y, mean_z
  real                          :: points(n,3)
!==============================================================================!

  ! First take aliasesd
  grid => mult % pnt_grid

  var_comb = HUGE
  var_comb_0 = 0.0
  vof_tol1 = MICRO
  vof_tol2 = MILI
  dist_cent = grid % vol(c) ** ONE_THIRD
  vof_int = HUGE

  n_tot = 0
  n_int = 0

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

  i_cell = .false.

  do while ( (n_tot < n) .and. (abs(var_comb - var_comb_0) > vof_tol1 .or. &
             (var_comb / dist_cent) > vof_tol2) )
    i_cell = .false.
    ! check if p is inside cell
    do while (i_cell .eqv. .false.)
      call random_number(r_num)
      p(1) = xmin + (xmax-xmin) * r_num
      call random_number(r_num)
      p(2) = ymin+ (ymax-ymin) * r_num
      call random_number(r_num)
      p(3) = zmin+ (zmax-zmin) * r_num
      i_cell = Check_Inside_Cell(mult, c, p)
    end do
    n_tot = n_tot + 1

    points(n_tot,:) = p(:)

    ! check if p is inside function:

    res_func = sqrt( ((p(1)-cent_x)/radius_1)**2    &
                   + ((p(2)-cent_y)/radius_2)**2    &
                   + ((p(3)-cent_z)/radius_3)**2 )

    if (res_func <= 1.0) then
      n_int = n_int + 1
    end if

    ! Compute maximum difference c/r center:
    if (mod(n_tot,500) == 0 .and. n_tot > 0) then
      var_comb_0 = var_comb
      vof_int = real(n_int) / real(n_tot)
      mean_x = sum(points(1: n_tot,1)) / real(n_tot)
      mean_y = sum(points(1: n_tot,2)) / real(n_tot)
      mean_z = sum(points(1: n_tot,3)) / real(n_tot)
      var_comb = abs(mean_x-grid % xc(c))   &
               + abs(mean_y-grid % yc(c))   &
               + abs(mean_z-grid % zc(c))
    end if
  end do

  end subroutine
