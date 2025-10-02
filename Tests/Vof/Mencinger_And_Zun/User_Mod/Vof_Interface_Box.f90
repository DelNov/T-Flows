!==============================================================================!
  subroutine Vof_Interface_Box(Vof,      &
                               c,        &
                               n_xyz,    &
                               dd,       &
                               vof_int)
!------------------------------------------------------------------------------!
!   Computes volume fraction of cell at interface                              !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Vof_Type),      target :: Vof
  integer, intent(in)         :: c
  real                        :: n_xyz(6,3)
  real                        :: dd(6)
  real                        :: vof_int
!------------------------------[Local parameters]------------------------------!
  integer, parameter :: N = 10000
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
  logical                  :: l_cel
  integer                  :: nod, n_int, n_tot, fu
  integer                  :: ee, n_cylinders, i_vari, n_vari
  real                     :: r_num, res_func
  real                     :: xmin, xmax, ymin, ymax, zmin, zmax
  real, allocatable        :: p(:,:)
  real                     :: vof_0, vof_tol1, vof_tol2
  real                     :: avg_x, avg_y, avg_z
  real                     :: var_comb, var_comb_0, dist_cent
  real                     :: mean_x, mean_y, mean_z
  real                     :: points(N,3)
!==============================================================================!

  ! First take aliasesd
  Grid => Vof % pnt_grid

  ! Initialize variables
  var_comb   = HUGE
  var_comb_0 = 0.0
  vof_tol1   = MICRO
  vof_tol2   = MILI
  dist_cent  = Grid % vol(c) ** ONE_THIRD
  vof_int    = HUGE

  n_tot = 0
  n_int = 0
  allocate(p(1,3))

  ! Find bounding box:
  xmin =  HUGE; ymin =  HUGE; zmin =  HUGE;
  xmax = -HUGE; ymax = -HUGE; zmax = -HUGE;

  do nod = 1, Grid % cells_n_nodes(c)
    xmin = min(xmin, Grid % xn(Grid % cells_n(nod,c)))
    ymin = min(ymin, Grid % yn(Grid % cells_n(nod,c)))
    zmin = min(zmin, Grid % zn(Grid % cells_n(nod,c)))
    xmax = max(xmax, Grid % xn(Grid % cells_n(nod,c)))
    ymax = max(ymax, Grid % yn(Grid % cells_n(nod,c)))
    zmax = max(zmax, Grid % zn(Grid % cells_n(nod,c)))
  end do

  l_cel = .false.

  do while ( (n_tot < N) .and. (abs(var_comb - var_comb_0) > vof_tol1 .or. &
             (var_comb / dist_cent) > vof_tol2) )

    l_cel = .false.

    ! Check if p is inside cell
    do while (l_cel .eqv. .false.)
      call random_number(r_num)
      p(1,1) = xmin + (xmax-xmin) * r_num
      call random_number(r_num)
      p(1,2) = ymin+ (ymax-ymin) * r_num
      call random_number(r_num)
      p(1,3) = zmin+ (zmax-zmin) * r_num
      l_cel = Grid % Is_Point_In_Cell(c, p(1,1), p(1,2), p(1,3))
    end do
    n_tot = n_tot + 1

    points(n_tot,:) = p(1,:)

    ! Check if p is inside function:
    if(Check_Inside_Box(Vof, p, dd, n_xyz) == 1) then
      n_int = n_int + 1
    end if

    ! Compute maximum difference c/r center:
    if(mod(n_tot, 500) == 0 .and. n_tot > 0) then
      var_comb_0 = var_comb
      vof_int = real(n_int) / real(n_tot)
      mean_x = sum(points(1: n_tot,1)) / real(n_tot)
      mean_y = sum(points(1: n_tot,2)) / real(n_tot)
      mean_z = sum(points(1: n_tot,3)) / real(n_tot)
      var_comb = abs(mean_x-Grid % xc(c))   &
               + abs(mean_y-Grid % yc(c))   &
               + abs(mean_z-Grid % zc(c))
    end if
  end do

  end subroutine
