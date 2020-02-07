!==============================================================================!
  subroutine Turb_Mod_Src_Kin_K_Eps_Zeta_F(turb, sol)
!------------------------------------------------------------------------------!
!   Computes the source terms in kin transport equation.                       !
!------------------------------------------------------------------------------!
!   In kinetic energy equation there are two source terms:                     !
!                                                                              !
!     /
!    |                                                                         !
!    | (density (p_kin - eps)) dV                                              !
!    |                                                                         !
!   /                                                                          !
!                                                                              !
!------------------------------------------------------------------------------!
  implicit none
!--------------------------------[Arguments]-----------------------------------!
  type(Turb_Type),   target :: turb
  type(Solver_Type), target :: sol
!---------------------------------[Calling]------------------------------------!
  real :: Roughness_Coefficient
  real :: Tau_Wall_Low_Re
  real :: Tau_Wall_Rough_Walls
  real :: Y_Plus_Low_Re
  real :: Y_Plus_Rough_Walls
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type),  pointer :: flow
  type(Grid_Type),   pointer :: grid
  type(Var_Type),    pointer :: u, v, w, t
  type(Var_Type),    pointer :: kin, eps, zeta, f, ut, vt, wt, t2
  type(Matrix_Type), pointer :: a
  real,              pointer :: b(:)
  integer                    :: c, c1, c2, s
  real                       :: u_tan, u_tau
  real                       :: lf, ebf, p_kin_int, p_kin_wf, l_rans_d, l_rans_v
  real                       :: l_rans, l_sgs, u_rans, u_sgs, kin_vis
  real                       :: z_o, alpha_d, alpha_v, l_sgs_d, l_sgs_v 
  real                       :: ut_log_law, vt_log_law, wt_log_law
  real                       :: nx, ny, nz, qx, qy, qz, g_buoy_wall
!==============================================================================!
!   Dimensions:                                                                !
!                                                                              !
!   production    p_kin    [m^2/s^3]   | rate-of-strain  shear     [1/s]       !
!   dissipation   eps % n  [m^2/s^3]   | turb. visc.     vis_t     [kg/(m*s)]  !
!   wall shear s. tau_wall [kg/(m*s^2)]| dyn visc.       viscosity [kg/(m*s)]  !
!   density       density  [kg/m^3]    | turb. kin en.   kin % n   [m^2/s^2]   !
!   cell volume   vol      [m^3]       | length          lf        [m]         !
!   left hand s.  A        [kg/s]      | right hand s.   b         [kg*m^2/s^3]!
!------------------------------------------------------------------------------!
!   p_kin = 2*vis_t / density S_ij S_ij                                        !
!   shear = sqrt(2 S_ij S_ij)                                                  !
!------------------------------------------------------------------------------!

  ! Take aliases
  flow => turb % pnt_flow
  grid => flow % pnt_grid
  call Field_Mod_Alias_Momentum   (flow, u, v, w)
  call Field_Mod_Alias_Energy     (flow, t)
  call Turb_Mod_Alias_K_Eps_Zeta_F(turb, kin, eps, zeta, f)
  call Turb_Mod_Alias_Heat_Fluxes (turb, ut, vt, wt)
  call Solver_Mod_Alias_System    (sol, a, b)
  call Turb_Mod_Alias_T2          (turb, t2)

  ! Production source:
  do c = 1, grid % n_cells
    turb % p_kin(c) = turb % vis_t(c) * flow % shear(c)**2
    b(c) = b(c) + turb % p_kin(c) * grid % vol(c)
  end do

  if(buoyancy) then
    do c = 1, grid % n_cells
      turb % g_buoy(c) = -flow % beta             &
                         * (grav_x * ut % n(c) +  &
                            grav_y * vt % n(c) +  &
                            grav_z * wt % n(c))   &
                          * flow % density(c)
      b(c) = b(c) + max(0.0, turb % g_buoy(c) * grid % vol(c))
             a % val(a % dia(c)) = a % val(a % dia(c))         &
                                 + max(0.0,-turb % g_buoy(c)   &
                                 * grid % vol(c)               &
                                 / (kin % n(c) + TINY))
    end do
  end if

  if(turbulence_model .eq. HYBRID_LES_RANS) then
    do c = 1, grid % n_cells

      lf = grid % vol(c)**ONE_THIRD

      ! Distance switch
      l_sgs_d  = 0.8 * lf
      l_rans_d = 0.41 * grid % wall_dist(c)
      alpha_d  = max(1.0,l_rans_d/l_sgs_d)

      ! Velocity switch
      l_sgs_v  = lf * flow % shear(c)
      l_rans_v = sqrt(kin % n(c) * zeta % n(c))
      alpha_v  = l_rans_v/l_sgs_v

      if( (hybrid_les_rans_switch .eq. SWITCH_DISTANCE)  &
          .and. (alpha_d < 1.05)                         &
          .or.                                           &
          (hybrid_les_rans_switch .eq. SWITCH_VELOCITY)  &
          .and. (alpha_v < 0.85) ) then
        a % val(a % dia(c)) = a % val(a % dia(c))             &
                            + flow % density(c) * eps % n(c)  &
                            / (kin % n(c) + TINY) * grid % vol(c)
      else
        a % val(a % dia(c)) = a % val(a % dia(c))   &
          + flow % density(c)                       &
          * min(alpha_d**1.4 * eps % n(c), kin % n(c)**1.5 / (lf*0.01))  &
          / (kin % n(c) + TINY) * grid % vol(c)
      end if
    end do
  else  ! turbuence model will be K_EPS_ZETA_F
    do c = 1, grid % n_cells
      a % val(a % dia(c)) = a % val(a % dia(c))             &
                          + flow % density(c) * eps % n(c)  &
                          / (kin % n(c) + TINY) * grid % vol(c)

    end do
  end if

  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    if(c2 < 0) then
      ! Kinematic viscosities
      kin_vis = flow % viscosity(c1) / flow % density(c1)

      if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALL .or. &
         Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALLFL) then

        ! Compute tangential velocity component
        u_tan = Field_Mod_U_Tan(flow, s)

        u_tau = c_mu25 * sqrt(kin % n(c1))

        turb % y_plus(c1) = Y_Plus_Low_Re(turb,                  &
                                          u_tau,                 &
                                          grid % wall_dist(c1),  &
                                          kin_vis)

        turb % tau_wall(c1) = Tau_Wall_Low_Re(turb,               &
                                              flow % density(c1), &
                                              u_tau,              &
                                              u_tan,              &
                                              turb % y_plus(c1))

        ebf = Turb_Mod_Ebf_Momentum(turb, c1)

        p_kin_wf  = turb % tau_wall(c1) * c_mu25 * sqrt(kin % n(c1))  &
                  / (grid % wall_dist(c1) * kappa)

        p_kin_int = turb % vis_t(c1) * flow % shear(c1)**2

        turb % p_kin(c1) = p_kin_wf

        if(turb % rough_walls) then
          z_o = Roughness_Coefficient(turb, turb % z_o_f(c1))
          turb % y_plus(c1) = Y_Plus_Rough_Walls(turb,                  &
                                                 u_tau,                 &
                                                 grid % wall_dist(c1),  &
                                                 kin_vis)

          turb % tau_wall(c1) = Tau_Wall_Rough_Walls(turb,                  &
                                                     flow % density(c1),    &
                                                     u_tau,                 &
                                                     u_tan,                 &
                                                     grid % wall_dist(c1),  &
                                                     z_o)

          turb % p_kin(c1) = turb % tau_wall(c1) * c_mu25 * sqrt(kin % n(c1)) &
                           / (kappa*(grid % wall_dist(c1)+z_o))

        end if ! rough_walls

        b(c1) = b(c1) + (turb % p_kin(c1)  &
              - turb % vis_t(c1) * flow % shear(c1)**2) * grid % vol(c1)

        ! Implementation of wall function for buoyancy-driven flows
        if(buoyancy) then

          nx = grid % sx(s) / grid % s(s)
          ny = grid % sy(s) / grid % s(s)
          nz = grid % sz(s) / grid % s(s)
          qx = t % q(c2) * nx
          qy = t % q(c2) * ny
          qz = t % q(c2) * nz

          ut_log_law = - turb % con_w(c1) &
                     * (t % n(c2) - t % n(c1))/grid % wall_dist(c1) * nx
          vt_log_law = - turb % con_w(c1) &
                     * (t % n(c2) - t % n(c1))/grid % wall_dist(c1) * ny
          wt_log_law = - turb % con_w(c1) &
                     * (t % n(c2) - t % n(c1))/grid % wall_dist(c1) * nz

          ut % n(c1) = ut %n(c1)  * exp(-1.0 * ebf) &
                     + ut_log_law * exp(-1.0 / ebf)
          vt % n(c1) = vt %n(c1)  * exp(-1.0 * ebf) &
                     + vt_log_law * exp(-1.0 / ebf)
          wt % n(c1) = wt %n(c1)  * exp(-1.0 * ebf) &
                     + wt_log_law * exp(-1.0 / ebf)

          if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALL) &
          t % q(c2) = turb % con_w(c1) * (t % n(c1) &
                      - t % n(c2))/grid % wall_dist(c1)

          g_buoy_wall = flow % beta*abs(grav_x + grav_y + grav_z)  &
                        * sqrt(abs(t % q(c2))                      &
                        * c_mu_theta5                              &
                        * sqrt(abs(t2 % n(c1) * kin % n(c1))))

          ! This limiter is preventing unphysical solutions
          ! when flow starts to develop
          g_buoy_wall = min(3.0*turb % p_kin(c1),g_buoy_wall)

          ! Clean up b(c) from old values of g_buoy
          b(c1)      = b(c1) - turb % g_buoy(c1) * grid % vol(c1)

          turb % g_buoy(c1) = turb % g_buoy(c1) * exp(-1.0 * ebf) &
                     + g_buoy_wall * exp(-1.0 / ebf)

          ! Add new values of g_buoy based on wall function approach
          b(c1)      = b(c1) + turb % g_buoy(c1) * grid % vol(c1)

        end if ! buoyancy

      end if  ! Grid_Mod_Bnd_Cond_Type(grid,c2).eq.WALL or WALLFL
    end if    ! c2 < 0
  end do

  end subroutine
