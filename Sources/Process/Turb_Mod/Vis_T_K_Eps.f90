!==============================================================================!
  subroutine Turb_Mod_Vis_T_K_Eps(turb)
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
  type(Turb_Type), target :: turb
!---------------------------------[Calling]------------------------------------!
  real :: Roughness_Coefficient
  real :: Tau_Wall_Low_Re
  real :: U_Plus_Log_Law
  real :: U_Plus_Rough_Walls
  real :: Y_Plus_Low_Re
  real :: Y_Plus_Rough_Walls
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type), pointer :: flow
  type(Grid_Type),  pointer :: grid
  type(Var_Type),   pointer :: u, v, w
  type(Var_Type),   pointer :: kin, eps
  integer                   :: c1, c2, s, c
  real                      :: pr, beta, ebf
  real                      :: u_tan, u_tau
  real                      :: kin_vis, u_plus, y_star, re_t, f_mu
  real                      :: z_o
!==============================================================================!
!   Dimensions:                                                                !
!                                                                              !
!   production    p_kin    [m^2/s^3]   | rate-of-strain  shear    [1/s]        !
!   dissipation   eps % n  [m^2/s^3]   | turb. visc.     vis_t    [kg/(m*s)]   !
!   wall shear s. tau_wall [kg/(m*s^2)]| dyn visc.       viscosity[kg/(m*s)]   !
!   density       density  [kg/m^3]    | turb. kin en.   kin % n  [m^2/s^2]    !
!   cell volume   vol      [m^3]       | length          lf       [m]          !
!   left hand s.  A        [kg/s]      | right hand s.   b        [kg*m^2/s^3] !
!   wall visc.    vis_w    [kg/(m*s)]  | kinematic viscosity      [m^2/s]      !
!   thermal cap.  capacity[m^2/(s^2*K)]| therm. conductivity     [kg*m/(s^3*K)]!
!------------------------------------------------------------------------------!
!   p_kin = 2*vis_t / density S_ij S_ij                                        !
!   shear = sqrt(2 S_ij S_ij)                                                  !
!------------------------------------------------------------------------------!

  ! Take aliases
  flow => turb % pnt_flow
  grid => flow % pnt_grid
  call Field_Mod_Alias_Momentum(flow, u, v, w)
  call Turb_Mod_Alias_K_Eps    (turb, kin, eps)

  do c = 1, grid % n_cells

    ! Kinematic viscosities
    kin_vis = flow % viscosity(c) / flow % density(c)

    re_t =  flow % density(c) * kin % n(c)**2  &
         / (flow % viscosity(c) * eps % n(c))

    y_star = (kin_vis * eps % n(c))**0.25 * grid % wall_dist(c)/kin_vis

    f_mu = (1.0 -     exp(-y_star/14.0))**2   &
         * (1.0 + 5.0*exp(-(re_t/200.0) * (re_t/200.0) ) /re_t**0.75)

    f_mu = min(1.0,f_mu)

    turb % vis_t(c) = f_mu * c_mu * flow % density(c) * kin % n(c)**2  &
                    / eps % n(c)
  end do

  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    if(c2 < 0) then
      if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALL .or.  &
         Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALLFL) then

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

        u_plus = U_Plus_Log_Law(turb, turb % y_plus(c1))

        if(turb % y_plus(c1) < 3.0) then
          turb % vis_w(c1) = turb % vis_t(c1) + flow % viscosity(c1)
        else
          turb % vis_w(c1) =    turb % y_plus(c1) * flow % viscosity(c1)  &
                           / (  turb % y_plus(c1) * exp(-1.0 * ebf)      &
                              + u_plus * exp(-1.0/ebf) + TINY)
        end if

        turb % y_plus(c1) = Y_Plus_Low_Re(turb,                  &
                                          u_tau,                 &
                                          grid % wall_dist(c1),  &
                                          kin_vis)

        if(turb % rough_walls) then
          z_o = Roughness_Coefficient(turb, turb % z_o_f(c1))
          turb % y_plus(c1) = Y_Plus_Rough_Walls(turb,                  &
                                                 u_tau,                 &
                                                 grid % wall_dist(c1),  &
                                                 kin_vis)
          u_plus     = U_Plus_Rough_Walls(turb, grid % wall_dist(c1))
          turb % vis_w(c1) = turb % y_plus(c1) * flow % viscosity(c1) / u_plus
        end if

        if(heat_transfer) then
          pr   = Field_Mod_Prandtl_Number(flow, c1)
          pr_t = Turb_Mod_Prandtl_Number(turb, c1)
          beta = 9.24 * ((pr/pr_t)**0.75 - 1.0)  &
               * (1.0 + 0.28 * exp(-0.007*pr/pr_t))
          ebf = Turb_Mod_Ebf_Scalar(turb, c1, pr)
          turb % con_w(c1) =    turb % y_plus(c1)                         &
                              * flow % viscosity(c1)                      &
                              * flow % capacity(c1)                       &
                      / (  turb % y_plus(c1) * pr * exp(-1.0 * ebf)       &
                         + (u_plus + beta) * pr_t * exp(-1.0 / ebf) + TINY)
        end if
      end if  ! Grid_Mod_Bnd_Cond_Type(grid,c2).eq.WALL or WALLFL
    end if    ! c2 < 0
  end do

  call Grid_Mod_Exchange_Real(grid, turb % vis_t)
  call Grid_Mod_Exchange_Real(grid, turb % vis_w)
  if(heat_transfer) then
    call Grid_Mod_Exchange_Real(grid, turb % con_w)
  end if

  end subroutine
