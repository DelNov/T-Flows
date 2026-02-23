!==============================================================================!
  subroutine Src_Kin_K_Omega_Sst(Turb, Sol)
!------------------------------------------------------------------------------!
!   Computes the source terms in kin transport equation for k-omega-sst model  !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Turb_Type),  target :: Turb
  type(Solver_Type), target :: Sol
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type),  pointer :: Flow
  type(Grid_Type),   pointer :: Grid
  type(Var_Type),    pointer :: u, v, w, t
  type(Var_Type),    pointer :: kin, omega, ut, vt, wt, t2
  type(Matrix_Type), pointer :: A
  real,              pointer :: b(:)
  integer                    :: c, c1, c2, s, reg
  real                       :: u_tan, u_tau
  real                       :: kin_vis
  real                       :: p_kin_int, p_kin_wf
  real                       :: z_o, ebf
  real                       :: ut_log_law, vt_log_law, wt_log_law
  real                       :: nx, ny, nz, qx, qy, qz, g_buoy_wall
!------------------------------------------------------------------------------!
!   Dimensions:                                                                !
!                                                                              !
!   production    p_kin    [m^2/s^3]   | rate-of-strain  shear     [1/s]       !
!   dissipation   eps % n  [m^2/s^3]   | turb. visc.     vis_t     [kg/(m*s)]  !
!   wall shear s. tau_wall [kg/(m*s^2)]| dyn visc.       viscosity [kg/(m*s)]  !
!   density       density  [kg/m^3]    | turb. kin en.   kin % n   [m^2/s^2]   !
!   cell volume   vol      [m^3]       | length          lf        [m]         !
!   left hand s.  A        [kg/s]      | kin. visc       kin_vis   [m^2/s]     !
!   right hand s. b        [kg*m^2/s^3]|                                       !
!------------------------------------------------------------------------------!
!   p_kin = 2*vis_t / density S_ij S_ij                                        !
!   shear = sqrt(2 S_ij S_ij)                                                  !
!==============================================================================!

  ! Take aliases
  Flow => Turb % pnt_flow
  Grid => Flow % pnt_grid
  call Flow % Alias_Momentum   (u, v, w)
  call Flow % Alias_Energy     (t)
  call Turb % Alias_Heat_Fluxes(ut, vt, wt)
  call Turb % Alias_T2         (t2)
  call Sol % Alias_Native      (A, b)

  kin   => turb % kin
  omega => turb % omega
  !-----------------------------------------!
  !   Compute the sources in the interior   !
  !-----------------------------------------!

  ! Production source:
  do c = Cells_In_Domain()
    Turb % p_kin(c) = Turb % vis_t(c) * Flow % shear(c)**2
    b(c) = b(c) + Turb % p_kin(c) * Grid % vol(c)

    A % val(A % dia(c)) = A % val(A % dia(c)) + &
         Flow % density(c) * Turb % beta_star * omega % n(c) * Grid % vol(c)

    if(Flow % buoyancy .eq. THERMALLY_DRIVEN) then
      Turb % g_buoy(c) = -Flow % beta                   &
                       * (  Flow % grav_x * ut % n(c)   &
                          + Flow % grav_y * vt % n(c)   &
                          + Flow % grav_z * wt % n(c))  &
                       * Flow % density(c)
      b(c) = b(c) + max(0.0, Turb % g_buoy(c) * Grid % vol(c))
      A % val(A % dia(c)) = A % val(A % dia(c))        &
                          + max(0.0,-Turb % g_buoy(c)  &
                          * Grid % vol(c) / (kin % n(c) + TINY))
    end if
  end do

  !-----------------------------------------------!
  !  Compute the sources in the near wall cells   !
  !-----------------------------------------------!
  do reg = Boundary_Regions()
    if(Grid % region % type(reg) .eq. WALL .or.  &
       Grid % region % type(reg) .eq. WALLFL) then
      do s = Faces_In_Region(reg)
        c1 = Grid % faces_c(1,s)
        c2 = Grid % faces_c(2,s)

        Assert(c2 < 0)  ! just to be on the safe side

        ! Kinematic viscosity
        kin_vis = Flow % viscosity(c1) / Flow % density(c1)

        ! Set up roughness coefficient
        z_o = Turb % Roughness_Coeff(c1, c2)

        ! Compute tangential velocity component
        u_tan = Flow % U_Tan(s)

        u_tau = Turb % c_mu25 * sqrt(kin % n(c1))

        Turb % y_plus(c1) = Turb % Y_Plus_Rough_Walls(    &
                                   u_tau,                 &
                                   Grid % wall_dist(c1),  &
                                   kin_vis,               &
                                   z_o)

        Turb % tau_wall(c1) = Turb % Tau_Wall_Log_Law(              &
                                              Flow % density(c1),   &
                                              u_tau,                &
                                              u_tan,                &
                                              Grid % wall_dist(c1), &
                                              Turb % y_plus(c1),    &
                                              z_o)

        ebf = Turb % Ebf_Momentum(c1)

        p_kin_wf  = Turb % tau_wall(c1) * Turb % c_mu25 * sqrt(kin % n(c1))  &
                  / ((Grid % wall_dist(c1) + z_o) * Turb % kappa)

        p_kin_int = Turb % vis_t(c1) * Flow % shear(c1)**2

        Turb % p_kin(c1) = exp(-1.0 * ebf) * p_kin_int   &
                         + exp(-1.0 / ebf) * p_kin_wf

        b(c1) = b(c1)                                                        &
              + (Turb % p_kin(c1) - Turb % vis_t(c1) * Flow % shear(c1)**2)  &
              * Grid % vol(c1)

        ! Implementation of wall function for buoyancy-driven flows
        if(Flow % buoyancy .eq. THERMALLY_DRIVEN) then

          nx = Grid % sx(s) / Grid % s(s)
          ny = Grid % sy(s) / Grid % s(s)
          nz = Grid % sz(s) / Grid % s(s)
          qx = t % q(c2) * nx
          qy = t % q(c2) * ny
          qz = t % q(c2) * nz

          ut_log_law = - Turb % con_w(c1)                                 &
                     / (Flow % density(c1) * Flow % capacity(c1))         &
                     * (t % n(c2) - t % n(c1))/Grid % wall_dist(c1) * nx
          vt_log_law = - Turb % con_w(c1)                                 &
                     / (Flow % density(c1) * Flow % capacity(c1))         &
                     * (t % n(c2) - t % n(c1))/Grid % wall_dist(c1) * ny
          wt_log_law = - Turb % con_w(c1)                                 &
                     / (Flow % density(c1) * Flow % capacity(c1))         &
                     * (t % n(c2) - t % n(c1))/Grid % wall_dist(c1) * nz

          ut % n(c1) = ut % n(c1) * exp(-1.0 * ebf)  &
                     + ut_log_law * exp(-1.0 / ebf)
          vt % n(c1) = vt % n(c1) * exp(-1.0 * ebf)  &
                     + vt_log_law * exp(-1.0 / ebf)
          wt % n(c1) = wt % n(c1) * exp(-1.0 * ebf)  &
                     + wt_log_law * exp(-1.0 / ebf)

          if(Grid % Bnd_Cond_Type(c2) .eq. WALL)                    &
            t % q(c2) = Turb % con_w(c1) * (t % n(c1) - t % n(c2))  &
                      / Grid % wall_dist(c1)

          g_buoy_wall = Flow % density(c1)                                    &
                      * Flow % beta * abs(  Flow % grav_x                     &
                                          + Flow % grav_y                     &
                                          + Flow % grav_z)                    &
                      * sqrt(abs(  t % q(c2)                                  &
                                 / (Flow % density(c1)*Flow % capacity(c1)))  &
                      * Turb % c_mu_theta5                                    &
                      * sqrt(abs(t2 % n(c1) * kin % n(c1))))

          ! Clean up b(c) from old values of g_buoy
          b(c1)      = b(c1) - Turb % g_buoy(c1) * Grid % vol(c1)

          Turb % g_buoy(c1) = Turb % g_buoy(c1) * exp(-1.0 * ebf) &
                            + g_buoy_wall * exp(-1.0 / ebf)

          ! Add new values of g_buoy based on wall function approach
          b(c1)      = b(c1) + Turb % g_buoy(c1) * Grid % vol(c1)

        end if  ! Flow % buoyancy .eq. THERMALLY_DRIVEN
      end do    ! faces in regions
    end if      ! region is WALL or WALLFL
  end do        ! through regions

  call Grid % Exchange_Cells_Real(kin % n)

  end subroutine
