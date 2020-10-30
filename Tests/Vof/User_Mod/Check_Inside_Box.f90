!==============================================================================!
  integer function Check_Inside_Box(mult, p, dd, n_xyz, c)
!------------------------------------------------------------------------------!
!   Determine if node of cell c lies inside box. Alteratively, if c is not     !
!   given it can be used to determine if any given point p is inside the box   !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Multiphase_Type), target :: mult
  real            , allocatable :: p(:,:)
  real                          :: dd(6)
  real                          :: n_xyz(6,3)
  integer            , optional :: c
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),      pointer :: grid
  integer                       :: i_fac,n
  integer                       :: infac_box(6)
  integer, allocatable          :: innod_box(:)
  integer                       :: sum_inside
  real                          :: res_dummy, dist, norm
!==============================================================================!

  ! First take aliasesd
  grid => mult % pnt_grid

  ! loop in cell faces:
  Check_Inside_Box = 0

  if (present(c)) then

    allocate(innod_box(grid % cells_n_nodes(c)))

    innod_box = 0

    ! for every node:
    do n = 1, grid % cells_n_nodes(c)
      infac_box = 0
      do i_fac = 1, 6
        res_dummy = n_xyz(i_fac,1) * grid % xn(grid % cells_n(n,c))      &
                  + n_xyz(i_fac,2) * grid % yn(grid % cells_n(n,c))      &
                  + n_xyz(i_fac,3) * grid % zn(grid % cells_n(n,c))
        norm = norm2(n_xyz(i_fac,1:3))
        dist = dot_product(n_xyz(i_fac,1:3)                &
             / norm, (/( grid % xn(grid % cells_n(n,c))    &
                       - p(i_fac,1) ),                     &
                       ( grid % yn(grid % cells_n(n,c))    &
                       - p(i_fac,2) ),                     &
                       ( grid % zn(grid % cells_n(n,c))    &
                       - p(i_fac,3) )/) )
        if (res_dummy <= dd(i_fac)) then
          infac_box(i_fac) = 1
        end if
      end do
      if (sum(infac_box) == 6) then
        innod_box(n) = 1
      end if
    end do

    sum_inside = sum(innod_box)
    if (sum_inside == grid % cells_n_nodes(c)) then
      Check_Inside_Box = 2
    else if (sum_inside < grid % cells_n_nodes(c) .and.   &
             sum_inside > 0) then
      Check_Inside_Box = 1
    end if
  else
    infac_box = 0
    do i_fac =1, 6
      res_dummy = n_xyz(i_fac,1) * p(1,1)     &
                + n_xyz(i_fac,2) * p(1,2)     &
                + n_xyz(i_fac,3) * p(1,3)
      if (res_dummy <= dd(i_fac)) then
        infac_box(i_fac) = 1
      end if
    end do

    if (sum(infac_box) == 6) then
      Check_Inside_Box = 1
    end if
  end if

  end function
