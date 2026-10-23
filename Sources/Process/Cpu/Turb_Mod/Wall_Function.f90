!==============================================================================!
  subroutine Wall_Function(Turb)
!------------------------------------------------------------------------------!
!   Computes turbulent viscosity for near-wall regions with wall functions     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Turb_Type), target :: Turb
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type), pointer :: Flow
  type(Grid_Type),  pointer :: Grid
  type(Var_Type),   pointer :: u, v, w, t, kin
  type(Var_Type),   pointer :: vis
  integer                   :: s, c1, c2, reg
  real                      :: u_tau, u_tan, nu
  real                      :: beta, pr, ebf, u_plus, sc, z_o, kin_vis
  real                      :: pr_tf, sc_tf
  integer                   :: it
  real                      :: a_mol, g_beta, d_t, t_c, y_x, th_x, th_log
  real                      :: g_nor, g_tan, con_nat, xi, w_nat
!------------------------------[Local parameters]------------------------------!
  real, parameter           :: A_POW = 8.3
  real, parameter           :: B_POW = 1.0/7.0
  real, parameter           :: C_NAT = 0.427  ! Hoelling and Herwig, J. Fluid
  real, parameter           :: D_NAT = 1.93   ! Mech., 541, pp. 383-397, 2005
  real, parameter           :: XI_LO = 0.1    ! below it, log-law only
  real, parameter           :: XI_HI = 0.2    ! above it, buoyant law only
!==============================================================================!

  ! Take aliases
  Flow => Turb % pnt_flow
  Grid => Flow % pnt_grid
  vis  => Turb % vis
  kin  => Turb % kin
  call Flow % Alias_Momentum(u, v, w)
  call Flow % Alias_Energy  (t)

  !-------------------!
  !   Wall function   !
  !-------------------+--------------!
  !   Law of the wall:               !
  !                                  !
  !   u+ = yz+  for z+ < 11.81       !
  !                                  !
  !   and                            !
  !                                  !
  !   u+ = A(y+)^B   for y+ > 11.81  !
  !                                  !
  !   with: A = 8.3 and B = 1/7      !
  !                                  !
  !----------------------------------+----------!
  !   The procedure below should be activated   !
  !   only if wall function approach is used.   !
  !---------------------------------------------!
  do reg = Boundary_Regions()

    if(Grid % region % type(reg) .eq. WALL .or.  &
       Grid % region % type(reg) .eq. WALLFL) then

      do s = Faces_In_Region(reg)
        c1 = Grid % faces_c(1,s)
        c2 = Grid % faces_c(2,s)

        kin_vis =  Flow % viscosity(c1) / Flow % density(c1)

        ! Set up roughness coefficient
        z_o = Turb % Roughness_Coeff(c1, c2)

        u_tan = Flow % U_Tan(s)

        nu = Flow % viscosity(c1) / Flow % density(c1)

        ! Calculate u_tau for smooth wall
        if(Turb % model == K_EPS .or. Turb % model == K_EPS_ZETA_F .or. &
           Turb % model == HYBRID_LES_RANS .or. Turb % model == K_OMEGA_SST) then
          u_tau = Turb % c_mu25 * sqrt(kin % n(c1))
        else
          u_tau = Turb % U_Tau_Log_Law(u_tan,                &
                                       Grid % wall_dist(c1), &
                                       kin_vis,              &
                                       z_o)
        end if

        ! Calculate u_tau according to Monin-Obukov Similarity Theory
        if(Flow % heat_transfer .and. Turb % monin_obukov) then
          u_tau = u_tau * Turb % Monin_Obukov_Momentum(abs(u_tan),  &
                  Grid % wall_dist(c1), z_o, t % n(c1),             &
                  t % n(c2), abs(Flow % grav_z))
        end if

        ! Calculate y+
        Turb % y_plus(c1) = Turb % Y_Plus_Rough_Walls(    &
                                   u_tau,                 &
                                   Grid % wall_dist(c1),  &
                                   kin_vis,               &
                                   z_o)

        u_plus = Turb % U_Plus_Log_Law(               &
                               Grid % wall_dist(c1),  &
                               Turb % y_plus(c1),     &
                               z_o)


        ebf = Turb % Ebf_Momentum(c1)

        if(Turb % y_plus(c1) < 3.0) then
          Turb % vis_w(c1) = Flow % viscosity(c1)
        else
          Turb % vis_w(c1) =    Turb % y_plus(c1) * Flow % viscosity(c1)  &
                           / (  Turb % y_plus(c1) * exp(-1.0 * ebf)       &
                              + u_plus * exp(-1.0/ebf) + TINY)
        end if

        if(Flow % heat_transfer) then
          ! Turbulent Prandtl number.  For LES models Prandtl_Turb gives
          ! the one for sub-grid scale fluxes, but the wall function
          ! bridges the layer below the first cell centre, where the
          ! log law holds, so pr_t from the control file is used for them
          pr_tf = Turb % Prandtl_Turb(c1)
          if(Turb % Les()) pr_tf = pr_t
          pr   = Flow % Prandtl_Numb(c1)          ! laminar Prandtl number
          beta = Turb % Beta_Scalar(pr, pr_tf)
          ! According to Toparlar et al. 2019 paper
          ! "CFD simulation of the near-neutral atmospheric boundary layer:
          ! New temperature inlet profile consistent with wall functions"
          if(z_o .gt. TINY) then
            beta = 0.0
          end if

          ! For y+ tending to zero, the expression below tends to molecular
          ! conductivity, but gives zero for y+ exactly zero.  With fluid at
          ! rest (y+ is zero at all walls) no heat could then enter through
          ! a wall with prescribed temperature, and a flow driven by buoyancy
          ! alone would never start.  Hence, the limit is set explicitly.
          if(Turb % y_plus(c1) < NANO) then
            Turb % con_w(c1) = Flow % conductivity(c1)
          else
            ebf = Turb % Ebf_Scalar(c1, pr)
            Turb % con_w(c1) =    Turb % y_plus(c1)                         &
                                * Flow % viscosity(c1)                      &
                                * Flow % capacity(c1)                       &
                        / (  Turb % y_plus(c1) * pr * exp(-1.0 * ebf)       &
                           + (u_plus + beta) * pr_tf * exp(-1.0 / ebf) + TINY)
          end if

          !--------------------------------------------------------------!
          !   Temperature wall function for natural convection, scaled   !
          !   with buoyancy instead of friction velocity, according to:  !
          !                                                              !
          !     Hoelling, M. and Herwig, H., "Asymptotic analysis of     !
          !     the near-wall region of turbulent natural convection     !
          !     flows", J. Fluid Mech., vol. 541, pp. 383-397, 2005.     !
          !                                                              !
          !     theta_x = y_x                    in viscous sublayer     !
          !     theta_x = C ln(y_x) + D          in overlap layer        !
          !                                                              !
          !   with: theta_x = (t_w - t) / t_c                            !
          !         y_x     = y * (g * beta * t_c / a^2)^(1/3)           !
          !                                                              !
          !   The two are blended with exponent 3, and t_c follows       !
          !   from temperature difference between wall and cell centre   !
          !                                                              !
          !   It is used only where buoyancy dominates over shear.  The  !
          !   criterion is the ratio of the two wall length scales:      !
          !                                                              !
          !     xi = y_x / y+ = (nu / u_tau) / (a^2/(g beta t_c))^(1/3)  !
          !                                                              !
          !   whose cube is Richardson number based on wall quantities.  !
          !   Below XI_LO standard wall function is kept, above XI_HI    !
          !   the buoyant one is used, with linear blend in between.     !
          !                                                              !
          !   At the walls with stable stratification, only the          !
          !   component of gravity tangential to the wall is used.       !
          !--------------------------------------------------------------!
          if((Turb % model .eq. HYBRID_LES_RANS .or.   &
              Turb % model .eq. K_EPS_ZETA_F)   .and.  &
             Flow % buoyancy .eq. THERMALLY_DRIVEN) then
            a_mol  = Flow % conductivity(c1)  &
                   / (Flow % density(c1) * Flow % capacity(c1))
            d_t    = abs(t % n(c2) - t % n(c1))

            ! Gravity; normal (positive out of the domain) and tangential
            g_beta = sqrt(  Flow % grav_x**2  &
                          + Flow % grav_y**2  &
                          + Flow % grav_z**2)
            g_nor  = (  Flow % grav_x * Grid % sx(s)  &
                      + Flow % grav_y * Grid % sy(s)  &
                      + Flow % grav_z * Grid % sz(s)) / Grid % s(s)
            g_tan  = sqrt(max(g_beta**2 - g_nor**2, 0.0))

            ! Stable stratification: warm wall above or cold wall below fluid
            if(g_nor * (t % n(c2) - t % n(c1)) * Flow % beta < 0.0) then
              g_beta = g_tan
            end if
            g_beta = abs(Flow % beta) * g_beta

            !----------------------------------------------!
            !   Natural convection part (buoyant scaling)  !
            !----------------------------------------------!
            con_nat = Flow % conductivity(c1)

            if(d_t > TINY .and. g_beta > TINY) then

              ! Initial guess from the linear profile
              t_c = (a_mol**2 * (d_t / Grid % wall_dist(c1))**3 / g_beta)**0.25

              do it = 1, 10
                y_x    = Grid % wall_dist(c1)  &
                       * (g_beta * t_c / a_mol**2)**ONE_THIRD
                th_x   = y_x
                th_log = C_NAT * log(y_x) + D_NAT
                if(th_log > 0.0) then
                  th_x = (1.0 / y_x**3 + 1.0 / th_log**3)**(-ONE_THIRD)
                end if
                t_c = d_t / th_x
              end do

              con_nat = Flow % conductivity(c1) * y_x / th_x
            end if

            !-----------------------------------------------------!
            !   Blend with the standard wall function (in con_w)  !
            !-----------------------------------------------------!
            xi = 0.0
            if(d_t > TINY .and. g_beta > TINY) then
              xi = y_x / max(Turb % y_plus(c1), TINY)
            end if
            w_nat = (xi - XI_LO) / (XI_HI - XI_LO)
            w_nat = max(0.0, min(1.0, w_nat))

            Turb % con_w(c1) = (1.0 - w_nat) * Turb % con_w(c1)  &
                             +        w_nat  * con_nat
            Turb % w_nat(c1) = w_nat
          end if

          if(Turb % monin_obukov) then
            Turb % con_w(c1) = pr_tf * Turb % con_w(c1)                 &
                             * Turb % Monin_Obukov_Thermal(abs(u_tan),  &
                               Grid % wall_dist(c1), z_o, t % n(c1),    &
                               t % n(c2),                               &
                               abs(Flow % grav_z))
          end if
        end if  ! heat transfer

        if(Flow % n_scalars > 0) then
          sc   = Flow % Schmidt_Numb(c1)          ! laminar Schmidt number

          ! Turbulent Schmidt number from the control file, also for LES
          ! models, for the same reason as pr_tf above
          sc_tf = sc_t

          beta = Turb % Beta_Scalar(sc, sc_tf)
          ! According to Toparlar et al. 2019 paper
          ! "CFD simulation of the near-neutral atmospheric boundary layer:
          ! New temperature inlet profile consistent with wall functions"
          if(z_o .gt. TINY) then
            beta = 0.0
          end if

          ! The same limit as for con_w above: molecular diffusivity
          if(Turb % y_plus(c1) < NANO) then
            Turb % diff_w(c1) = (Flow % viscosity(c1)/Flow % density(c1)) / sc
          else
            ebf = Turb % Ebf_Scalar(c1, sc)
            Turb % diff_w(c1) =  Turb % y_plus(c1)                  &
                * (Flow % viscosity(c1)/Flow % density(c1))         &
                / (Turb % y_plus(c1) * sc * exp(-1.0 * ebf)         &
                 + (u_plus + beta) * sc_tf * exp(-1.0 / ebf) + TINY)
          end if

          if(Turb % monin_obukov) then
            Turb % diff_w(c1) = sc_tf * Turb % diff_w(c1)                     &
                           * Turb % Monin_Obukov_Thermal(abs(u_tan),          &
                             Grid % wall_dist(c1), z_o, t % n(c1), t % n(c2), &
                             abs(Flow % grav_z))
          end if
        end if  ! n_scalars > 0

      end do
    end if      ! Grid % region % type(reg) .eq. WALL
  end do

  call Grid % Exchange_Cells_Real(Turb % y_plus)
  call Grid % Exchange_Cells_Real(Turb % w_nat)
  call Grid % Exchange_Cells_Real(Turb % vis_w)
  if(Flow % n_scalars > 0) then
    call Grid % Exchange_Cells_Real(Turb % diff_w)
  end if
  if(Flow % heat_transfer) then
    call Grid % Exchange_Cells_Real(Turb % con_w)
  end if

  end subroutine
