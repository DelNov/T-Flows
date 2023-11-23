!==============================================================================!
  logical function Is_Point_In_Cell(Grid, c, xp, yp, zp)
!------------------------------------------------------------------------------!
!>  This function checks whether a given point (specified by coordinates xp, yp,
!>  zp) lies within a specified cell (c) in a computational grid. It is a
!>  geometric function used to assess the spatial relationship between a point
!>  and the grid's cells.
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * Iterating over cell faces:                                               !
!     - The function examines each face of the cell, determining its normal    !
!       vector and the center point.                                           !
!   * Determining face orientation:                                            !
!     - It calculates the orientation of the face's normal vector based on the !
!       cell's position relative to its neighboring cells.                     !
!   * Checking point's position relative to faces:                             !
!     - For each face, the function computes a dot product between the face's  !
!       normal vector and the vector from the face's center to the point.      !
!   * Counting inward-facing normals:                                          !
!     - If the dot product is negative or zero, it implies the point lies      !
!       towards the inside of the face (or on it). The function counts such    !
!       instances.                                                             !
!   * Final determination:                                                     !
!     - If the count of inward-facing normals equals the number of cell faces, !
!       the function concludes that the point is inside the cell and returns   !
!       true, else false.                                                      !
!------------------------------------------------------------------------------!
  class(Grid_Type), target :: Grid        !! computational grid
  integer, intent(in)      :: c           !! cell index
  real,    intent(in)      :: xp, yp, zp  !! point's cordinates
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
