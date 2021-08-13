!==============================================================================!
!  Compile with:                                                               !
!                                                                              !
!  gfortran -o Sucking_Problem -fdefault-real-8 ./Sucking_Problem.f90          !
!------------------------------------------------------------------------------!

!==============================================================================!
  module Properties_Mod
!------------------------------------------------------------------------------!
!  Module holding all physical properties                                      !
!------------------------------------------------------------------------------!
    real, parameter :: RHO_G   =    0.597    ! [kg/m^3]
    real, parameter :: RHO_L   =  958.4      ! [kg/m^3]
    real, parameter :: K_G     =    0.025    ! [W/m/K]
    real, parameter :: K_L     =    0.679    ! [W/m/K]
    real, parameter :: CP_G    = 2030.0      ! [J/kg K]
    real, parameter :: CP_L    = 4216.0      ! [J/kg K]
    real, parameter :: H_LG    =    2.26e+6  ! [J/kg]
    real, parameter :: T_WALL  =   10.00     ! [K]
    real, parameter :: T_SAT   =   10.00     ! [K]
    real, parameter :: T_INF   =   15.00     ! [K]
    real, parameter :: ALPHA_G = K_G / (RHO_G * CP_G)
    real, parameter :: ALPHA_L = K_L / (RHO_L * CP_L)
    real, parameter :: PI = 4.0 * atan(1.0)
  end module

!==============================================================================!
  real function Transc(beta)
!------------------------------------------------------------------------------!
!   Equation (31) from Rajkotwala et al. (2019)                                !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Properties_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real :: beta
!-----------------------------------[Locals]-----------------------------------!
  real :: inb  ! in braces
!==============================================================================!

  ! Term in braces; appears twice
  inb = beta**2 * (ALPHA_G * RHO_G**2)  &
                / (ALPHA_L * RHO_L**2)

  Transc = (beta                                                          &
            - ( (T_INF-T_SAT) * CP_G * K_L * sqrt(ALPHA_G) * exp(-inb) )  &
            / (    H_LG * K_G * sqrt(PI * ALPHA_L) * erfc(sqrt(inb))   )  &
           )

  end function

!==============================================================================!
  real function Beta(beta_l, beta_r, tol)
!------------------------------------------------------------------------------!
!   Solution of equation (31) from Rajkotwala et al. (2019)                    !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Properties_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real :: beta_l, beta_r, tol
!----------------------------------[Calling]-----------------------------------!
  real :: Transc
!-----------------------------------[Locals]-----------------------------------!
  integer :: i
  real    :: xm_l, xm_r, rhs, rhs_l, rhs_r
!==============================================================================!

  print '(a)', '#==============================='
  print '(a)', '# Iterative calculation of Beta'
  print '(a)', '#-------------------------------'

  rhs = 0.0

  print '(a, es15.7)', '# Initial rhs = ', rhs
  print '(a)',         '#  Left value,    Right value,' //  &
                       '   Left residual, Right residual'

  ! Set initial values of moving boundary
  xm_l = beta_l
  xm_r = beta_r

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
  print '(a, es15.7)', '# Convergence reached, beta = ', 0.5 * (xm_l + xm_r)
  Beta = 0.5 * (xm_l + xm_r)

  end function

!==============================================================================!
  subroutine Temperature_Distribution(beta, n, tim, ipos, fpos)
!------------------------------------------------------------------------------!
!   Equation (52) from M. Irfan, M. Muradoglu (2017)                           !
!   Note that equation (32) from Rajkotwala et al. (2019) has an error         !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Properties_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real    :: beta  ! beta sent from the main function
  integer :: n     ! time step (needed to form file name)
  real    :: tim   ! current time
  real    :: ipos  ! interface position
  real    :: fpos  ! final position
!-----------------------------------[Locals]-----------------------------------!
  integer            :: i
  real               :: cpos, opos, tcur, told
  character(32)      :: fname = 'temperature_distribution_XXX.dat'
  character(29)      :: gname = 'temperature_gradients_XXX.dat'
  integer, parameter :: TD = 13  ! temperature distribution
  integer, parameter :: TG = 14  ! temperature gradients
!==============================================================================!

  write(fname(26:28), '(i3.3)') n
  write(gname(23:25), '(i3.3)') n

  print '(a,f12.5)', '# Time:          ', tim
  print '(a,a)', '# Creating file: ', fname
  print '(a,a)', '# Creating file: ', gname
  open(file=fname, unit=TD)
  open(file=gname, unit=TG)

  do i = 0, 1000

    ! Calculate new position
    cpos = 1.5 * fpos * real(i)/1000.0  ! cpos will range from 0 to 4.8 * fpos

    ! Store temperatures
    if(cpos > ipos) then
      tcur =  T_INF                                 &
           - (T_INF-T_SAT)                          &
           / erfc(beta * (RHO_G * sqrt(ALPHA_G))    &
                        /(RHO_L * sqrt(ALPHA_L))    &
                 )                                  &
           * erfc(  cpos / (2.0*sqrt(ALPHA_L*tim))   &
                  + beta * (RHO_G - RHO_L)/RHO_L * sqrt(ALPHA_G/ALPHA_L) )
    else
      tcur = T_SAT
    end if
    write(TD,*) cpos, tcur

    ! Store temperature gradients
    if(i > 0) then
      write(TG,*)  (cpos+opos)*0.5, (tcur-told)/(cpos-opos)
    end if

    ! Store previous position
    opos = cpos
    told = tcur
  end do

  close(TD)
  close(TG)

  end subroutine

!==============================================================================!
  program Sucking_Problem_Main
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Properties_Mod
!------------------------------------------------------------------------------!
  implicit none
!----------------------------------[Calling]-----------------------------------!
  real :: Beta
!-----------------------------------[Locals]-----------------------------------!
  integer            :: n, n_dt
  integer            :: i
  real               :: beta_main
  real               :: ipos, opos, start_pos, final_pos, dt
  real               :: start_time, final_time, time
  character(22)      :: fname = 'interface_position.dat'
  character(22)      :: gname = 'interface_velocity.dat'
  integer, parameter :: IP = 11  ! interface position
  integer, parameter :: IV = 12  ! interface velocity
!==============================================================================!

  print '(a)', '#============================='
  print '(a)', '#'
  print '(a)', '# Solving the Sucking problem'
  print '(a)', '#'
  print '(a)', '#-----------------------------'

  ! Calculate main xi (initial guesses and tolerance)
  beta_main = Beta(0.0,  1.0,  1.0e-10)

  ! Calculate position
  start_time =   0.1                              ! starting time        [s]
  final_time =   0.6                              ! final time           [s]
  n_dt       = 100                                ! number of time steps [1]
  dt         =  (final_time - start_time) / n_dt  ! time step            [s]

  ! Final interface position
  start_pos = 2.0 * beta_main * sqrt(ALPHA_G * start_time)
  final_pos = 2.0 * beta_main * sqrt(ALPHA_G * final_time)
  print '(a)', '#=============================================='
  print '(a,es15.7)', '# Start interface position is: ', start_pos
  print '(a,es15.7)', '# Final interface position is: ', final_pos
  print '(a)', '#----------------------------------------------'

  print '(a,a)', '# Creating file: ', fname
  print '(a,a)', '# Creating file: ', gname
  open(file=fname, unit=IP)
  open(file=gname, unit=IV)

  write(IP, '(a)')  '# Time [s]       Position [m]'
  write(IV, '(a)')  '# Time [s]       Velocity [m/s]'
  opos = 0.0                                     ! old position
  do n = 0, n_dt

    ! Advance time and position
    time = start_time + dt * n                     ! current time
    ipos = 2.0 * beta_main * sqrt(ALPHA_G * time)  ! current position

    ! Save temperature distribution
    if(mod(n,10) .eq. 0) then
      call Temperature_Distribution(beta_main, n, time, ipos, final_pos)
    end if

    write(IP, '(2es15.7)') time, ipos
    if(n > 0) then
      write(IV, '(2es15.7)') time - dt*0.5, (ipos-opos)/dt
    end if
    opos = ipos
  end do

  close(IP)
  close(IV)

  end program
