!==============================================================================!
  subroutine Src_T2(Turb, Sol)
!------------------------------------------------------------------------------!
!   Computes the source terms in t2 transport equation for k-eps_t2 model      !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Turb_Type),  target :: Turb
  type(Solver_Type), target :: Sol
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type),  pointer :: Flow
  type(Grid_Type),   pointer :: Grid
  type(Var_Type),    pointer :: u, v, w, t
  type(Var_Type),    pointer :: kin, eps, ut, vt, wt, t2, omega
  type(Matrix_Type), pointer :: A
  real,              pointer :: b(:)
  integer                    :: c, c1, c2, s, reg
  real                       :: kin_vis, p_t2_wall, ebf, u_tau
  real                       :: ut_sgdh, vt_sgdh, wt_sgdh, z_o, eps_eff
!------------------------------------------------------------------------------!
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
!==============================================================================!

  ! Take aliases
  Flow => Turb % pnt_flow
  Grid => Flow % pnt_grid
  call Flow % Alias_Momentum   (u, v, w)
  call Flow % Alias_Energy     (t)
  call Turb % Alias_K_Eps      (kin, eps)
  call Turb % Alias_Heat_Fluxes(ut, vt, wt)
  call Turb % Alias_T2         (t2)
  call Sol % Alias_Native      (A, b)
  omega => Turb % omega

  !-----------------------------------------!
  !   Compute the sources in all the cells  !
  !-----------------------------------------!

  ! Temperature gradients are needed
  call Flow % Grad_Variable(t)

  ! Production source:
  do c = Cells_In_Domain_And_Buffers()

    !-------------------------------------------------------------------!
    !   ut, vt and wt defined by AFM or GGDH could lead to divergence   !
    !-------------------------------------------------------------------!
    pr_t = max(Turb % Prandtl_Turb(c), TINY)
    ut_sgdh = -Turb % vis_t(c) / Flow % density(c) / pr_t * t % x(c)
    vt_sgdh = -Turb % vis_t(c) / Flow % density(c) / pr_t * t % y(c)
    wt_sgdh = -Turb % vis_t(c) / Flow % density(c) / pr_t * t % z(c)

    if(Turb % model .eq. HYBRID_LES_RANS) then
      ut_sgdh = -Turb % vis_t_eff(c) / Flow % density(c) / pr_t * t % x(c)
      vt_sgdh = -Turb % vis_t_eff(c) / Flow % density(c) / pr_t * t % y(c)
      wt_sgdh = -Turb % vis_t_eff(c) / Flow % density(c) / pr_t * t % z(c)
    end if

    Turb % p_t2(c) = - 2.0 * Flow % density(c)       &
                           * (  ut_sgdh * t % x(c)   &
                           +    vt_sgdh * t % y(c)   &
                           +    wt_sgdh * t % z(c))

    b(c) = b(c) + Turb % p_t2(c) * Grid % vol(c)

  ! Negative contribution

    if (Turb % model == K_OMEGA_SST) then
      eps_eff = Turb % beta_star * kin % n(c) * omega % n(c)
    else
      eps_eff = eps % n(c)
    end if

    A % val(A % dia(c)) = A % val(A % dia(c)) +  &
         2.0 * Flow % density(c) * eps_eff  &
             / (kin % n(c) + TINY) * Grid % vol(c)

  end do

  !------------------------------------------------------------------------!
  !   Implementation of wall function approach for buoyancy-driven flows   !
  !------------------------------------------------------------------------!
  do reg = Boundary_Regions()
    if(Grid % region % type(reg) .eq. WALL .or.  &
       Grid % region % type(reg) .eq. WALLFL) then
      do s = Faces_In_Region(reg)
        c1 = Grid % faces_c(1,s)
        c2 = Grid % faces_c(2,s)

        ! Set up roughness coefficient
        z_o = Turb % Roughness_Coeff(c1, c2)

        ! Kinematic viscosities
        kin_vis = Flow % viscosity(c1) / Flow % density(c1)

        u_tau = Turb % c_mu25 * sqrt(kin % n(c1))

        Turb % y_plus(c1) = Turb % Y_Plus_Rough_Walls(u_tau,                &
                                                      Grid % wall_dist(c1), &
                                                      kin_vis,              &
                                                      z_o)
        ebf = Turb % Ebf_Momentum(c1)

        p_t2_wall = Flow % density(c1)                                       &
                  * abs(t % q(c2)/(Flow % density(c1)*Flow % capacity(c1)))  &
                  * Turb % c_mu_theta5*sqrt(abs(t2 % n(c1)))                 &
                  / (Turb % kappa_theta * Turb % c_mu25 * Grid % wall_dist(c1))

        b(c1) = b(c1) - Turb % p_t2(c1) * Grid % vol(c1)

        if(Turb % y_plus(c1) > 11.0) then
          Turb % p_t2(c1) = p_t2_wall
        else
          Turb % p_t2(c1) = (  Turb % p_t2(c1) * exp(-1.0 * ebf)  &
                                   + p_t2_wall * exp(-1.0 / ebf))
        end if

        b(c1) = b(c1) + Turb % p_t2(c1) * Grid % vol(c1)

        t2 % n(c2) = 0.0

      end do  ! faces in regions
    end if    ! region is WALL or WALLFL
  end do      ! through regions

  end subroutine

