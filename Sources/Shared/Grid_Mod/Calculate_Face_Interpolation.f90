!==============================================================================!
  subroutine Grid_Mod_Calculate_Face_Interpolation(grid)
!------------------------------------------------------------------------------!
!   Calculate interpolation factors for faces.                                 !
!   Should not be called from ProcessR, but from Genearate and Convert         !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer :: s, c1, c2, cnt
  real    :: xc1, yc1, zc1, xc2, yc2, zc2, nx, ny, nz, s_tot
  real    :: xi, yi, zi, lx, ly, lz, dsc1, dsc2
  real    :: f_avg, f_max, f_min, r_avg, r_max, r_min, r
!==============================================================================!

  !----------------!
  !   Old method   !
  !----------------!
  ! do s = 1, grid % n_faces
  !   c1 = grid % faces_c(1,s)
  !   c2 = grid % faces_c(2,s)
  !
  !   ! First cell
  !   xc1  = grid % xc(c1)
  !   yc1  = grid % yc(c1)
  !   zc1  = grid % zc(c1)
  !   dsc1 = Math_Mod_Distance(xc1, yc1, zc1,   &
  !                            grid % xf(s), grid % yf(s), grid % zf(s))
  !
  !   ! Second cell (pls. check if xsi=xc on the boundary)
  !   xc2  = grid % xc(c2) + grid % dx(s)
  !   yc2  = grid % yc(c2) + grid % dy(s)
  !   zc2  = grid % zc(c2) + grid % dz(s)
  !   dsc2 = Math_Mod_Distance(xc2, yc2, zc2,   &
  !                            grid % xf(s), grid % yf(s), grid % zf(s))
  !
  !   ! Interpolation factor
  !   grid % f(s) = dsc2 / (dsc1 + dsc2)
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

  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    ! Surface normal
    s_tot = sqrt(grid % sx(s)**2 + grid % sy(s)**2 + grid % sz(s)**2)
    nx = grid % sx(s) / s_tot
    ny = grid % sy(s) / s_tot
    nz = grid % sz(s) / s_tot

    ! First cell
    xc1  = grid % xc(c1)
    yc1  = grid % yc(c1)
    zc1  = grid % zc(c1)

    ! Second cell (with correction for periodicity)
    xc2  = grid % xc(c2) + grid % dx(s)
    yc2  = grid % yc(c2) + grid % dy(s)
    zc2  = grid % zc(c2) + grid % dz(s)

    ! Vector connecting cell centers c1 and c2
    lx = xc2 - xc1
    ly = yc2 - yc1
    lz = zc2 - zc1

    ! Distance from c1 to intersection
    dsc1 = (  (grid % xf(s) - xc1) * nx     &
            + (grid % yf(s) - yc1) * ny     &
            + (grid % zf(s) - zc1) * nz  )  &
         / (lx * nx + ly * ny + lz * nz)

    ! Intersection point
    xi = xc1 + dsc1 * lx
    yi = yc1 + dsc1 * ly
    zi = zc1 + dsc1 * lz

    grid % xr(s) = grid % xf(s) - xi
    grid % yr(s) = grid % yf(s) - yi
    grid % zr(s) = grid % zf(s) - zi

    ! Interpolation factor
    ! (Equation 2.19 in Denner's thesis)
    dsc2 = Math_Mod_Distance(xc2, yc2, zc2, xi, yi, zi)

    grid % f(s) = dsc2 / (dsc1 + dsc2)

    if(c2 > 0) then
      cnt = cnt + 1
      f_max = max(f_max, grid % f(s))
      f_min = min(f_min, grid % f(s))
      f_avg = f_avg + grid % f(s)
      r = sqrt(grid % xr(s)**2 + grid % yr(s)**2 + grid % zr(s)**2)  &
        / (dsc1 + dsc2)
      r_max = max(r_max, r)
      r_min = min(r_min, r)
      r_avg = r_avg + r
    end if
  end do

  !-----------------------------------!
  !   Print some info on the screen   !
  !-----------------------------------!
  f_avg = f_avg / cnt
  r_avg = r_avg / cnt
  print *, '#================================================================'
  print *, '# Face interpolation factors have been calculated'
  print *, '# Factors closer to 0.5 and shifts closer to 0.0 are better'
  print *, '#----------------------------------------------------------------'
  print '(a19,es12.5)', ' # Minimum factor: ', f_min
  print '(a19,es12.5)', ' # Maximum factor: ', f_max
  print '(a19,es12.5)', ' # Average factor: ', f_avg
  print '(a19,es12.5)', ' # Minimum shift : ', r_min
  print '(a19,es12.5)', ' # Maximum shift : ', r_max
  print '(a19,es12.5)', ' # Average shift : ', r_avg

  end subroutine
