!==============================================================================!
  program Parabolic
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
  real    :: r    = 0.51
  real    :: bulk = 1.0/1.25
  integer :: i
  integer :: n = 51

  real    :: r_l, r_r, r_c, y_l, y_r, y_c, dr, integral, area
!==============================================================================!

  print *, '#=================='
  print *, '# Number of points '
  print *, '#=================='
  print *, n

  print *, '#================================='
  print *, '#   Coordinate      Velocity'
  print *, '#================================='

  ! Distance between two nodes
  dr = r / (n-1)

  ! Print profile and compute integral
  area     = 0.0
  integral = 0.0
  do i = 1, n
    r_l = dr * (i-1)  ! r left
    r_r = dr *  i     ! r rigth
    r_c = 0.5 * (r_l + r_r)
    y_l = (r**2 - r_l**2) * bulk
    y_r = (r**2 - r_r**2) * bulk

    print '(2es16.5)', r-r_l, y_l

    if(i < n) then
      y_c      = (y_l + y_r) * 0.5
      area     = area     + dr * r_c
      integral = integral + y_c * dr * r_c
    end if
  end do

  print *, '# Bulk velocity is: ', integral / area

  ! Compute integral

  end program
