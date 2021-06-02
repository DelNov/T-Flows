!==============================================================================!
  logical function Is_Point_In_Cell(Grid, c, xp, yp, zp)
!------------------------------------------------------------------------------!
!   Determine is point x_p, y_p, z_p lies inside cell c                        !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type), target :: Grid
  integer, intent(in)      :: c
  real,    intent(in)      :: xp, yp, zp
!-----------------------------------[Locals]-----------------------------------!
  integer :: c1, c2, s, i_fac, cnt_in
  real    :: nx, ny, nz, xf, yf, zf, dot_prod
!==============================================================================!

  ! loop in cell faces:
  Is_Point_In_Cell = .false.

  cnt_in = 0

  do i_fac = 1, Grid % cells_n_faces(c)
    s  = Grid % cells_f(i_fac, c)
    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)

    ! Normal is positive from c1 to c2, so if c is c1, this is how it should be
    nx = Grid % sx(s) / Grid % s(s)
    ny = Grid % sy(s) / Grid % s(s)
    nz = Grid % sz(s) / Grid % s(s)

    xf = Grid % xf(s)
    yf = Grid % yf(s)
    zf = Grid % zf(s)

    ! Yet, of is c is actually c2, normal changes the sign (clear enough)
    if(c == c2) then
      nx = -nx
      ny = -ny
      nz = -nz
      if(Grid % faces_s(s) > 0) then
        xf = Grid % xf(Grid % faces_s(s))
        yf = Grid % yf(Grid % faces_s(s))
        zf = Grid % zf(Grid % faces_s(s))
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

  if(cnt_in == Grid % cells_n_faces(c)) then
    Is_Point_In_Cell = .true.
  end if

  end function
