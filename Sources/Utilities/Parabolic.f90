!==============================================================================!
  program Parabolic
!------------------------------------------------------------------------------!
!   Creates a parabolic velocity profile useful for laminar inflows            !
!                                                                              !
!   Compile with: gfortran -o Parabolic Parabolic.f90                          !
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
  real          :: x_s, x_e, bulk
  integer       :: n, i
  real          :: x_l, x_r, x_c, y_l, y_r, y_c, dx, lx, integral, area
  character(80) :: arg
!==============================================================================!

  ! Proper invocation
  if(command_argument_count() .eq. 4) then

    ! Get x_s, x_e and bulk
    call get_command_argument(1, arg);  read(arg, *) x_s
    call get_command_argument(2, arg);  read(arg, *) x_e
    call get_command_argument(3, arg);  read(arg, *) bulk
    call get_command_argument(4, arg);  read(arg, *) n

  ! Wrong invocation
  else
    print '(a)', '# You failed to invoke the program properly.'
    print '(a)', '# Correct invocation:'
    print '(a)', './Parabolic  x_start  x_end  bulk_velocity  n_points'
    stop
  end if

  print '(a)', '#    Number of points:'
  print '(i6)', n

  print '(a)', '#    Coordinate:     Velocity:'

  ! distance between two nodes
  lx = (x_e - x_s)
  dx = lx / (n-1)

  ! Print profile and compute integral
  area     = 0.0
  integral = 0.0
  do i = 1, n
    ! x_l and x_r are in range -0.5 -> +0.5
    x_l  = -0.5 + real((i-1))/real(n-1)  ! x left
    x_r  = -0.5 + real( i   )/real(n-1)  ! x rigth
    y_l = (1.0 - x_l ** 2 / 0.25) * bulk *  1.5
    y_r = (1.0 - x_r ** 2 / 0.25) * bulk *  1.5

    print '(2es16.5)', x_s + (i-1)*dx,  y_l

    if(i < n) then
      y_c      = (y_l + y_r) * 0.5
      area     = area     + dx
      integral = integral + y_c * dx
    end if
  end do

  ! Print integral
  print '(a,es13.3)', '# Bulk velocity: ', integral / area

  end program
