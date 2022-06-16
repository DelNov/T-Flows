!==============================================================================!
  subroutine Vis_T_K_Eps(Turb)
!------------------------------------------------------------------------------!
!   Computes the turbulent viscosity for RANS models.                          !
!                                                                              !
!   In the domain:                                                             !
!   For k-eps model :                                                          !
!                                                                              !
!   vis_t = c_mu * rho * k^2 * eps                                             !
!                                                                              !
!   On the boundary (wall viscosity):                                          !
!   vis_tw = y^+ * vis_t kappa / (E * ln(y^+))                                 !
!                                                                              !
!   For k-eps-v2f model :                                                      !
!                                                                              !
!   vis_t = CmuD * rho * Tsc  * vv                                             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Turb_Type), target :: Turb
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type), pointer :: Flow
  type(Grid_Type),  pointer :: Grid
  type(Var_Type),   pointer :: u, v, w
  type(Var_Type),   pointer :: kin, eps
  integer                   :: c1, c2, s, c
  real                      :: pr, beta, ebf, sc
  real                      :: u_tan, u_tau
  real                      :: kin_vis, u_plus, y_star, re_t, f_mu
  real                      :: z_o
!==============================================================================!
!   Dimensions:                                                                !
!                                                                              !
!   production    p_kin    [m^2/s^3]   | rate-of-strain  shear    [1/s]        !
!   dissipation   eps % n  [m^2/s^3]   | Turb. visc.     vis_t    [kg/(m*s)]   !
!   wall shear s. tau_wall [kg/(m*s^2)]| dyn visc.       viscosity[kg/(m*s)]   !
!   density       density  [kg/m^3]    | Turb. kin en.   kin % n  [m^2/s^2]    !
!   cell volume   vol      [m^3]       | length          lf       [m]          !
!   left hand s.  A        [kg/s]      | right hand s.   b        [kg*m^2/s^3] !
!   wall visc.    vis_w    [kg/(m*s)]  | kinematic viscosity      [m^2/s]      !
!   thermal cap.  capacity[m^2/(s^2*K)]| therm. conductivity     [kg*m/(s^3*K)]!
!------------------------------------------------------------------------------!
!   p_kin = 2*vis_t / density S_ij S_ij                                        !
!   shear = sqrt(2 S_ij S_ij)                                                  !
!------------------------------------------------------------------------------!

  ! Take aliases
  Flow => Turb % pnt_flow
  Grid => Flow % pnt_grid
  call Flow % Alias_Momentum(u, v, w)
  call Turb % Alias_K_Eps   (kin, eps)

  do c = 1, Grid % n_cells

    ! Kinematic viscosities
    kin_vis = Flow % viscosity(c) / Flow % density(c)

    re_t =  Flow % density(c) * kin % n(c)**2  &
         / (Flow % viscosity(c) * eps % n(c))

    y_star = (kin_vis * eps % n(c))**0.25 * Grid % wall_dist(c)/kin_vis

    f_mu = (1.0 -     exp(-y_star/14.0))**2   &
         * (1.0 + 5.0*exp(-(re_t/200.0) * (re_t/200.0) ) /re_t**0.75)

    f_mu = min(1.0,f_mu)

    Turb % vis_t(c) = f_mu * c_mu * Flow % density(c) * kin % n(c)**2  &
                      / (eps % n(c) + TINY)
  end do

  do s = 1, Grid % n_faces
    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)

    if(c2 < 0) then
      if(Grid % Bnd_Cond_Type(c2) .eq. WALL .or.  &
         Grid % Bnd_Cond_Type(c2) .eq. WALLFL) then

        kin_vis =  Flow % viscosity(c1) / Flow % density(c1)

        ! Set up roughness coefficient
        z_o = Turb % Roughness_Coefficient(Turb % z_o_f(c1))
        if(Turb % rough_walls) then
          z_o = max(Grid % wall_dist(c1)  &
              / (e_log * max(Turb % y_plus(c1), 1.0)), z_o)
        end if

        ! Compute tangential velocity component
        u_tan = Flow % U_Tan(s)

        u_tau = c_mu25 * sqrt(kin % n(c1))

        Turb % y_plus(c1) = Turb % Y_Plus_Rough_Walls(   &
                                   u_tau,                &
                                   Grid % wall_dist(c1), &
                                   kin_vis,              &
                                   z_o)

        Turb % tau_wall(c1) = Turb % Tau_Wall_Log_Law(              &
                                              Flow % density(c1),   &
                                              u_tau,                &
                                              u_tan,                &
                                              Grid % wall_dist(c1), &
                                              Turb % y_plus(c1),    &
                                              z_o)

        ebf = Turb % Ebf_Momentum(c1)

        u_plus = Turb % U_Plus_Log_Law(               &
                               Grid % wall_dist(c1),  &
                               Turb % y_plus(c1),     &
                               z_o)

        if(Turb % y_plus(c1) < 3.0) then
          Turb % vis_w(c1) = Turb % vis_t(c1) + Flow % viscosity(c1)
        else

          if(Turb % y_plus(c1) < 11.3 ) then
            ebf = 0.00001
          end if

          Turb % vis_w(c1) =    Turb % y_plus(c1) * Flow % viscosity(c1)  &
                           / (  Turb % y_plus(c1) * exp(-1.0 * ebf)       &
                              + u_plus * exp(-1.0/ebf) + TINY)

        end if

        if(Flow % heat_transfer) then
          pr   = Flow % Prandtl_Numb(c1)
          pr_t = Turb % Prandtl_Turb(c1)
          beta = 9.24 * ((pr/pr_t)**0.75 - 1.0)  &
               * (1.0 + 0.28 * exp(-0.007*pr/pr_t))
          ebf = Turb % Ebf_Scalar(c1, pr)
          Turb % con_w(c1) =    Turb % y_plus(c1)                         &
                              * Flow % viscosity(c1)                      &
                              * Flow % capacity(c1)                       &
                      / (  Turb % y_plus(c1) * pr * exp(-1.0 * ebf)       &
                         + (u_plus + beta) * pr_t * exp(-1.0 / ebf) + TINY)
        end if

        if(Flow % n_scalars > 0) then
          sc   = Flow % Schmidt_Numb(c1)            ! laminar Schmidt number
          beta = 9.24 * ((sc/sc_t)**0.75 - 1.0)                     &
               * (1.0 + 0.28 * exp(-0.007*sc/sc_t))
          ebf  = 0.01 * (sc * Turb % y_plus(c1)**4                  &
               / ((1.0 + 5.0 * sc**3 * Turb % y_plus(c1)) + TINY))
          Turb % diff_w(c1) =  Turb % y_plus(c1)                    &
               * (Flow % viscosity(c1)/Flow % density(c1))          &
               / (  Turb % y_plus(c1) * sc * exp(-1.0 * ebf)        &
               + (u_plus + beta) * sc_t * exp(-1.0 / ebf) + TINY)
        end if

      end if  ! Grid % Bnd_Cond_Type(c2).eq.WALL or WALLFL
    end if    ! c2 < 0
  end do

  call Grid % Exchange_Cells_Real(Turb % vis_t)
  call Grid % Exchange_Cells_Real(Turb % vis_w)
  if(Flow % heat_transfer) then
    call Grid % Exchange_Cells_Real(Turb % con_w)
  end if
  if(Flow % n_scalars > 0) then
    call Grid % Exchange_Cells_Real(Turb % diff_w)
  end if

  end subroutine
