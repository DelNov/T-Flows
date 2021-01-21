!==============================================================================!
  logical function Grid_Mod_Is_Point_In_Cell(grid, c, xp, yp, zp)
!------------------------------------------------------------------------------!
!   Determine is point x_p, y_p, z_p lies inside cell c                        !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type), target :: grid
  integer                 :: c
  real                    :: xp, yp, zp
!-----------------------------------[Locals]-----------------------------------!
  integer :: c1, c2, s, i_fac, cnt_in
  integer :: inside_c(MAX_CELLS_N_FACES)
  real    :: nx, ny, nz, xf, yf, zf, dot_prod
!==============================================================================!

  ! loop in cell faces:
  Grid_Mod_Is_Point_In_Cell = .false.

  cnt_in = 0

  do i_fac = 1, grid % cells_n_faces(c)
    s  = grid % cells_f(i_fac, c)
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    ! Normal is positive from c1 to c2, so if c is c1, this is how it should be
    nx = grid % sx(s) / grid % s(s)
    ny = grid % sy(s) / grid % s(s)
    nz = grid % sz(s) / grid % s(s)

    xf = grid % xf(s)
    yf = grid % yf(s)
    zf = grid % zf(s)

    ! Yet, of is c is actually c2, normal changes the sign (clear enough)
    if(c == c2) then
      nx = -nx
      ny = -ny
      nz = -nz
      if(grid % faces_s(s) > 0) then
        xf = grid % xf(grid % faces_s(s))
        yf = grid % yf(grid % faces_s(s))
        zf = grid % zf(grid % faces_s(s))
      end if
    end if

    ! Vector p - f is oriented from the face towards the point
    ! If point is inside the cell, dot_prod will turn out negative
    dot_prod = nx * (xp - xf)  &
             + ny * (yp - yf)  &
             + nz * (zp - zf)

    ! Count number of times point qualifies to be inside of the cell c
    if(dot_prod <= 0.0) then
      cnt_in = cnt_in + 1
    end if

  end do

  if(cnt_in == grid % cells_n_faces(c)) then
    Grid_Mod_Is_Point_In_Cell = .true.
  end if

  end function
