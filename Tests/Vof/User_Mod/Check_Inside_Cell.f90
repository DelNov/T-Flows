!==============================================================================!
  logical function Check_Inside_Cell(mult, c, p)
!------------------------------------------------------------------------------!
!   Determine is point p lies inside cell c                                    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Multiphase_Type), target :: mult
  integer                       :: c
  real                          :: p(3)
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),      pointer :: grid
  integer                       :: c1, c2, s, i_fac
  integer,          allocatable :: inside_c(:)
  real                          :: n_unit(3), dist
  real                          :: corr_x, corr_y, corr_z
!==============================================================================!

  ! First take aliasesd
  grid => mult % pnt_grid

  ! loop in cell faces:
  Check_Inside_Cell = .false.

  allocate(inside_c(grid % cells_n_faces(c)))

  inside_c = 0

  do i_fac = 1, grid % cells_n_faces(c)
    s = grid % cells_f(i_fac, c)
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    corr_x = 0.0; corr_y = 0.0; corr_z = 0.0

    if(c == c1) then
      n_unit =-(/grid % sx(s) / grid % s(s),    &
                 grid % sy(s) / grid % s(s),    &
                 grid % sz(s) / grid % s(s)/)
    else
      n_unit = (/grid % sx(s) / grid % s(s),    &
                 grid % sy(s) / grid % s(s),    &
                 grid % sz(s) / grid % s(s)/)

      ! Correction for periodic faces:
      call Grid_Mod_Correction_Periodicity(grid, s, corr_x, corr_y, corr_z)

      ! This is to make sure the precision of points location inside cells
      if (abs(corr_x) < PICO) corr_x = 0.0
      if (abs(corr_y) < PICO) corr_y = 0.0
      if (abs(corr_z) < PICO) corr_z = 0.0
    end if

    dist = dot_product( n_unit,(/p(1) - (grid % xf(s) + corr_x),      &
                                 p(2) - (grid % yf(s) + corr_y),      &
                                 p(3) - (grid % zf(s) + corr_z)/) )

    if (dist >= 0.0) then
      inside_c(i_fac) = 1
    end if
  end do

  if (sum(inside_c) == grid % cells_n_faces(c)) then
    Check_Inside_Cell = .true.
  end if

  end function
