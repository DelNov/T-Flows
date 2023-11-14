!==============================================================================!
  subroutine Calculate_Face_Interpolation(Grid)
!------------------------------------------------------------------------------!
!>  Calculate interpolation factors for faces.
!------------------------------------------------------------------------------!
!   Note: It shouldn't be called from Process, but from Genearate and Convert  !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type) :: Grid
!-----------------------------------[Locals]-----------------------------------!
  integer :: s, c1, c2, cnt
  real    :: xc1, yc1, zc1, xc2, yc2, zc2, nx, ny, nz, s_tot
  real    :: xi, yi, zi, lx, ly, lz, dsc1, dsc2
  real    :: f_avg, f_max, f_min, r_avg, r_max, r_min, r
!==============================================================================!

  !----------------!
  !   Old method   !
  !----------------!
  ! do s = 1, Grid % n_faces
  !   c1 = Grid % faces_c(1,s)
  !   c2 = Grid % faces_c(2,s)
  !
  !   ! First cell
  !   xc1  = Grid % xc(c1)
  !   yc1  = Grid % yc(c1)
  !   zc1  = Grid % zc(c1)
  !   dsc1 = Math % Distance(xc1, yc1, zc1,   &
  !                            Grid % xf(s), Grid % yf(s), Grid % zf(s))
  !
  !   ! Second cell (pls. check if xsi=xc on the boundary)
  !   xc2  = Grid % xc(c2) + Grid % dx(s)
  !   yc2  = Grid % yc(c2) + Grid % dy(s)
  !   zc2  = Grid % zc(c2) + Grid % dz(s)
  !   dsc2 = Math % Distance(xc2, yc2, zc2,   &
  !                            Grid % xf(s), Grid % yf(s), Grid % zf(s))
  !
  !   ! Interpolation factor
  !   Grid % f(s) = dsc2 / (dsc1 + dsc2)
  ! end do

  !-----------------------------------------------------------!
  !   New method as proposed in Fabian Denner's thesis        !
  !                                                           !
  !   Line-face intersection follows procedure described in:  !
  !   https://en.wikipedia.org/wiki/Line-plane_intersection   !
  !-----------------------------------------------------------!
  f_max = -HUGE
  f_min = +HUGE
  f_avg =  0.0
  r_max = -HUGE
  r_min = +HUGE
  r_avg =  0.0
  cnt   =  0

  do s = 1, Grid % n_faces
    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)

    ! Surface normal
    s_tot = sqrt(Grid % sx(s)**2 + Grid % sy(s)**2 + Grid % sz(s)**2)
    nx = Grid % sx(s) / s_tot
    ny = Grid % sy(s) / s_tot
    nz = Grid % sz(s) / s_tot

    ! First cell
    xc1  = Grid % xc(c1)
    yc1  = Grid % yc(c1)
    zc1  = Grid % zc(c1)

    ! Second cell (with correction for periodicity)
    xc2  = Grid % xc(c2) + Grid % dx(s)
    yc2  = Grid % yc(c2) + Grid % dy(s)
    zc2  = Grid % zc(c2) + Grid % dz(s)

    ! Vector connecting cell centers c1 and c2
    lx = xc2 - xc1
    ly = yc2 - yc1
    lz = zc2 - zc1

    ! Distance from c1 to intersection
    dsc1 = (  (Grid % xf(s) - xc1) * nx     &
            + (Grid % yf(s) - yc1) * ny     &
            + (Grid % zf(s) - zc1) * nz  )  &
         / (lx * nx + ly * ny + lz * nz)

    ! Intersection point
    xi = xc1 + dsc1 * lx
    yi = yc1 + dsc1 * ly
    zi = zc1 + dsc1 * lz

    Grid % rx(s) = Grid % xf(s) - xi
    Grid % ry(s) = Grid % yf(s) - yi
    Grid % rz(s) = Grid % zf(s) - zi

    ! Interpolation factor
    ! (Equation 2.19 in Denner's thesis)
    dsc1 = Math % Distance(xc1, yc1, zc1, xi, yi, zi)
    dsc2 = Math % Distance(xc2, yc2, zc2, xi, yi, zi)

    Grid % f(s) = dsc2 / (dsc1 + dsc2)

    if(c2 > 0) then
      cnt = cnt + 1
      f_max = max(f_max, Grid % f(s))
      f_min = min(f_min, Grid % f(s))
      f_avg = f_avg + Grid % f(s)
      r = sqrt(Grid % rx(s)**2 + Grid % ry(s)**2 + Grid % rz(s)**2)  &
        / (dsc1 + dsc2)
      r_max = max(r_max, r)
      r_min = min(r_min, r)
      r_avg = r_avg + r
    end if
  end do

  !-----------------------------------!
  !   Print some info on the screen   !
  !-----------------------------------!
  f_avg = f_avg / real(cnt)
  r_avg = r_avg / real(cnt)
  call Message % Framed(54,                                         &
           'Face interpolation factors have been calculated',       &
           'It is useful to keep in mind that factors closer '  //  &
           'to 0.5 and shifts closer to 0.0 are better')
  print '(a19,es12.5)', ' # Minimum factor: ', f_min
  print '(a19,es12.5)', ' # Maximum factor: ', f_max
  print '(a19,es12.5)', ' # Average factor: ', f_avg
  print '(a19,es12.5)', ' # Minimum shift : ', r_min
  print '(a19,es12.5)', ' # Maximum shift : ', r_max
  print '(a19,es12.5)', ' # Average shift : ', r_avg

  end subroutine
