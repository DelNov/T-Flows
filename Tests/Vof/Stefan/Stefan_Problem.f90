!==============================================================================!
!  Compile with:                                                               !
!                                                                              !
!  gfortran -o Stefan_Problem -fdefault-real-8 ./Stefan_Problem.f90            !
!------------------------------------------------------------------------------!

!==============================================================================!
  module Properties_Mod
!------------------------------------------------------------------------------!
!  Module holding all physical properties                                      !
!------------------------------------------------------------------------------!
    real, parameter :: CP_V    = 2030.0      ! [J/kg K]
    real, parameter :: LAMDA_V =    0.025    ! [W/m K]
    real, parameter :: RHO_V   =    0.597    ! [W/m K]
    real, parameter :: L       =    2.26e+6  ! [J/kg]
    real, parameter :: T_WALL  =  110.00     ! [K]
    real, parameter :: T_SAT   =  100.00     ! [K]
    real, parameter :: ALPHA_V = LAMDA_V / (RHO_V * CP_V)
    real, parameter :: PI = 4.0 * atan(1.0)
  end module

!==============================================================================!
  real function Transc(xi)
!------------------------------------------------------------------------------!
!   Left-hand side of equation (29) from Sato and Niceno, 2013                 !
!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real :: xi
!==============================================================================!

  Transc = xi * exp((xi**2)) * erf(xi)

  end function

!==============================================================================!
  real function Xi(xi_l, xi_r, tol)
!------------------------------------------------------------------------------!
!   Iterative solution of equation (29) from Sato and Niceno, 2013             !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Properties_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real :: xi_l, xi_r, tol
!----------------------------------[Calling]-----------------------------------!
  real :: Transc
!-----------------------------------[Locals]-----------------------------------!
  integer :: i
  real    :: xm_l, xm_r, rhs, rhs_l, rhs_r
!==============================================================================!

  print '(a)', '#============================='
  print '(a)', '# Iterative calculation of Xi'
  print '(a)', '#-----------------------------'

  rhs = CP_V * (T_WALL - T_SAT) / (sqrt(PI) * L)

  print '(a, es15.7)', '# rhs = ', rhs

  ! Set initial values of moving boundary
  xm_l = xi_l
  xm_r = xi_r

  ! Start iterations
  do i = 1, 200

    rhs_l = Transc(xm_l)
    rhs_r = Transc(xm_r)

    ! Check if converged
    if(abs(rhs - rhs_l) < tol .and. abs(rhs - rhs_r) < tol) goto 1

    print '(a, 4es15.7)', '#', xm_l, xm_r, rhs_l-rhs, rhs_r-rhs

    if(rhs_r > rhs .and. Transc((0.5*(xm_l+xm_r))) > rhs) then
      xm_r = 0.5 * (xm_l + xm_r)
    end if

    if(rhs_l < rhs .and. Transc((0.5*(xm_l+xm_r))) < rhs) then
      xm_l = 0.5 * (xm_l + xm_r)
    end if

  end do

  ! If here, it failed to reach convergence
  print '(a)', '# Failed to converge in function "Xi"'
  stop

  ! Convergence reached
1 continue
  print '(a, es15.7)', '# Convergence reached, xi = ', 0.5 * (xm_l + xm_r)
  Xi = 0.5 * (xm_l + xm_r)

  end function

!==============================================================================!
  subroutine Temperature_Distribution(xi, n, tim, ipos, fpos)
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Properties_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real    :: xi    ! xi sent from the main function
  integer :: n     ! time step (needed to form file name)
  real    :: tim   ! current time
  real    :: ipos  ! interface position
  real    :: fpos  ! final position
!-----------------------------------[Locals]-----------------------------------!
  integer       :: i
  real          :: cpos, opos, tcur, told
  character(32) :: fname = 'temperature_distribution_XXX.dat'
  character(29) :: gname = 'temperature_gradients_XXX.dat'
!==============================================================================!

  write(fname(26:28), '(i3.3)') n
  write(gname(23:25), '(i3.3)') n

  print '(a,a)', '# Creating file: ', fname
  print '(a,a)', '# Creating file: ', gname
  open(file=fname, unit=13)
  open(file=gname, unit=14)

  do i = 0, 100

    ! Calculate new position
    cpos = 1.2 * fpos * real(i)/100.0  ! cpos will range from 0 to 1.2 * fpos

    ! Store temperatures
    if(cpos < ipos) then
      tcur = T_WALL  &
           + ((T_SAT-T_WALL) / erf(xi)) * erf(cpos / (2.*sqrt(ALPHA_V*tim)))
    else
      tcur = T_SAT
    end if
    write(13,*)  cpos, tcur

    ! Store temperature gradients
    if(i > 0) then
      write(14,*)  (cpos+opos)*0.5, (tcur-told)/(cpos-opos)
    end if

    ! Store previous position
    opos = cpos
    told = tcur
  end do

  close(13)
  close(14)

  end subroutine

!==============================================================================!
  program Stefan_Problem_Main
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Properties_Mod
!------------------------------------------------------------------------------!
  implicit none
!----------------------------------[Calling]-----------------------------------!
  real :: Xi
!-----------------------------------[Locals]-----------------------------------!
  integer       :: n, n_dt
  integer       :: i
  real          :: xi_main
  real          :: ipos, opos, final_pos, dt, start_time, final_time, time
  character(22) :: fname = 'interface_position.dat'
  character(22) :: gname = 'interface_velocity.dat'
!==============================================================================!

  print '(a)', '#============================================'
  print '(a)', '#'
  print '(a)', '# Solving the Stefan problem for steam-water'
  print '(a)', '#'
  print '(a)', '#--------------------------------------------'

  ! Calculate main xi (initial guesses and tolerance)
  xi_main = Xi(0.0,  1.0,  1.0e-10)

  ! Calculate position
  start_time =   0.6766                           ! starting time        [s]
  final_time =   9.6766                           ! final time           [s]
  n_dt       = 100                                ! number of time steps [1]
  dt         =  (final_time - start_time) / n_dt  ! time step            [s]

  ! Final interface position
  final_pos = 2.0 * xi_main * sqrt(ALPHA_V * final_time)
  print '(a)', '#=============================================='
  print '(a,es15.7)', '# Final interface position is: ', final_pos
  print '(a)', '#----------------------------------------------'

  print '(a,a)', '# Creating file: ', fname
  print '(a,a)', '# Creating file: ', gname
  open(file=fname, unit=11)
  open(file=gname, unit=12)

  write(11, '(a)')  '# Time [s]       Position [m]'
  write(12, '(a)')  '# Time [s]       Velocity [m/s]'
  opos = 0.0                                     ! old position
  do n = 0, n_dt

    ! Advance time and position
    time = start_time + dt * n                   ! current time
    ipos = 2.0 * xi_main * sqrt(ALPHA_V * time)  ! current position

    ! Save temperature distribution
    if(mod(n,10) .eq. 0) then
      call Temperature_Distribution(xi_main, n, time, ipos, final_pos)
    end if

    write(11, '(2es15.7)') time - start_time, ipos
    if(n > 0) then
      write(12, '(2es15.7)') time - dt*0.5 - start_time, (ipos-opos)/dt
    end if
    opos = ipos
  end do

  close(11)
  close(12)

  end program
