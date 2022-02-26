include '../User_Mod/T_Sat.f90'

!==============================================================================!
  subroutine User_Mod_Beginning_Of_Iteration(Flow, Turb, Vof, Swarm, n, time)
!------------------------------------------------------------------------------!
!   This function is called at the beginning of time step.                     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),      target :: Flow
  type(Turb_Type),       target :: Turb
  type(Vof_Type),        target :: Vof
  type(Swarm_Type),      target :: Swarm
  integer                       :: n     ! time step
  real                          :: time  ! physical time
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: Grid
  type(Var_Type),   pointer :: u, v, w, t, scalar
  real, contiguous, pointer :: dens(:)
  integer                   :: dataIndex, c, s, c1, c2
  real                      :: dens_air(9), dens_h2o(9), temperat_int(9)
  real                      :: weight, a0, a1, a2, a3, a4, a5, b0, b1, b2, c0
  real                      :: aa, bb, cc, a_tmp, m_h2o, m_air, M, pv, t_int
  real                      :: t_film, k_film, d_film, h_d, m_evap
!==============================================================================!

  ! Take aliases
  Grid   => Flow % pnt_grid
  t      => Flow % t
  scalar => Flow % scalar(1)

  ! Temperature (and salt) dependent density
  temperat_int = (/10.0, 20.0, 30.0, 40.0, 50.0,  60.0,  70.0,   80.0, 90.0 /)
  dens_air = (/1.246, 1.204, 1.166, 1.127, 1.097, 1.067, 1.034, 1.0, 0.9734 /)
  dens_h2o = (/999.7, 998.3, 995.7, 992.3, 988.0, 983.0, 978.0, 972.0, 965.0/)

  !--------------------------------!
  !   Operations in upper domain   !
  !--------------------------------!
  if(Grid % name .eq. 'UPPER_DOM') then

    do c = 1, Grid % n_cells
      if(t % n(c) .gt. 90.0 .or. t % n(c) .lt. 10.0) then
         error stop "temperature out of density interpolation range [10 - 90C]"
      endif

      dataIndex = floor(t % n(c) /10) 
      if(t % n(c) .eq. 90.0) dataIndex = 8
      weight = (t % n(c) - temperat_int(dataIndex))                   &
             / (temperat_int(dataIndex+1) - temperat_int(dataIndex))

      Flow % density(c) = (1.0 - weight) * dens_h2o(dataIndex)      &
                        +        weight  * dens_h2o(dataIndex + 1)

      ! Salt influence taken from Millero and Huang 2009
      a_tmp = scalar % n(c) * 1000 ! convert to g/kg

      a0 =  8.197247E-01
      a1 = -3.779454E-03
      a2 =  6.821795E-05
      a3 = -8.009571E-07
      a4 =  6.158885E-09
      a5 = -2.001919E-11

      b0 = -5.808305E-03
      b1 =  5.354872E-05
      b2 = -4.714602E-07

      c0 =  5.249266E-04

      if (t % n(c) .ge. 0.0 .and. t % n(c) .le. 90.0) then
        aa = a0 + a1 * t % n(c) + a2 * t % n(c)**2 + &
             a3 * t % n(c)**3 + a4 * t %n(c)**4 + &
             a5 * t % n(c)**5

        bb = b0 + b1 * t % n(c) + b2 * t % n(c)**2
        cc = c0

        Flow % density(c) = Flow % density(c) + aa * scalar % n(c) &
                                              + bb * scalar % n(c)**1.5 &
                                              + cc * scalar % n(c)**2
      else
        error stop 'Temperature value out of range for ' //  &
                   'calculation of salt dependent density [0-90C]'
      end if

    end do

  end if

  !--------------------------------!
  !   Operations in lower domain   !
  !--------------------------------!
  if(Grid % name .eq. 'LOWER_DOM') then

    do c = 1, Grid % n_cells
      if(t % n(c) .gt. 90.0 .or. t % n(c) .lt. 10.0) then
         error stop 'Temperature out of density interpolation range [10 - 90C]'
      endif

      dataIndex = floor(t % n(c) /10) 
      if(t % n(c) .eq. 90.0) dataIndex = 8
      weight = (t % n(c) - temperat_int(dataIndex)) / &
               (temperat_int(dataIndex+1) - temperat_int(dataIndex))

      Flow % density(c)  = (1.0 - weight) * dens_air(dataIndex)      &
                         +        weight  * dens_air(dataIndex + 1)
    end do

    ! As long as no film domain exists
    m_h2o  = 18.0e-3
    m_air  = 28.0e-3
    t_film = 20.0
    h_d    = 2454e3 ! J/kg
    k_film = 0.6009 ! W/mK
    d_film = 1e-3  ! m

    do s = 1, Grid % n_faces
      c1 = Grid % faces_c(1,s)
      c2 = Grid % faces_c(2,s) !-> each face 2 cells, right cell
      if (c2 .lt. 0) then ! at boundary

        ! Check boundary by its name
        if (Var_Mod_Bnd_Cond_Name(scalar,c2) .eq. 'BOTTOM_WALL') then
          M = 1.0 / ( (1.0 - scalar % n(c1)) / m_air  &
                           + scalar % n(c1)  / m_h2o )
          pv = scalar % n(c1) * M/m_h2o * 1e5
          !print * , pv!
          call t_sat(t_int, pv)

          ! Underrelaxation of condensation absolutely necessary!!
          t_int = t_int * 0.1 + 0.9 * t_int_prev(c1)

          t_int_prev(c1) = t_int
          !t % n(c2) = t_int
          m_evap = (Flow % conductivity(c1) &
                  / Grid % wall_dist(c1) * (t % n(c1)-t_int)  &
                 + k_film / d_film * (t_film - t_int)) / h_d

          ! scalar % q(c2) = m_evap * Grid % sz(s)
        end if
      end if
    end do

    ! Positive for evaporation, negative for condensation
    if(this_proc < 2) then
      print * , 'm_evap =', m_evap * 3600, ' kg/m²h '
      print * , 't_int = ', t_int, 'Celsius'
    end if

  end if

  end subroutine
