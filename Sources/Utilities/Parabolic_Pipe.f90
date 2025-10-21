!==============================================================================!
  program Parabolic_Pipe
!------------------------------------------------------------------------------!
!   Creates a parabolic velocity profile useful for laminar inflows            !
!                                                                              !
!   Compile with: gfortran -o Parabolic_Pipe Parabolic_Pipe.f90                !
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
  real, parameter :: r_min = 0.0
  real            :: r_max, bulk
  integer         :: n, i
  real            :: r_pos, y_l, y_r, y_c, dr, r, r_c, integral, area
  character(80)   :: arg
!==============================================================================!

  ! Proper invocation
  if(command_argument_count() .eq. 3) then

    ! Get r_max, bulk and number of points
    call get_command_argument(1, arg);  read(arg, *) r_max
    call get_command_argument(2, arg);  read(arg, *) bulk
    call get_command_argument(3, arg);  read(arg, *) n

  ! Wrong invocation
  else
    print '(a)', '# You failed to invoke the program properly.'
    print '(a)', '# Correct invocation:'
    print '(a)', './Parabolic_Pipe  r_max  bulk_velocity  n_points'
    stop
  end if

  print '(a)', '#    Number of points:'
  print '(i6)', n

  print '(a)', '#    Coordinate:     Velocity:'

  ! Distance between two nodes
  r  = (r_max - r_min)
  dr = r / (n-1)

  ! Print profile and compute integral
  area     = 0.0
  integral = 0.0
  do i = n, 1, -1

    ! Current radial position
    r_pos = r_min + (i-1)*dr

    ! Pipe flow parabolic profile: u(r) = 2.0 * bulk * (1 - (r/r_max)^2)
    y_l = (1.0 - ( r_pos       / r_max)**2) * bulk * 2.0
    y_r = (1.0 - ((r_pos + dr) / r_max)**2) * bulk * 2.0

    print '(2es16.5)', r_pos, y_l

    if(i < n) then
      ! For integration, use next point
      y_c = (y_l + y_r) * 0.5
      r_c = r_pos + 0.5 * dr
      area     = area     +       r_c * dr
      integral = integral + y_c * r_c * dr
    end if
  end do

  ! Print integral
  print '(a,es13.3)', '# Bulk velocity: ', integral / area

  end program
