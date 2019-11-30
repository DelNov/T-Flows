!==============================================================================!
  subroutine Turb_Mod_Vis_T_K_Eps_Zeta_F(turb)
!------------------------------------------------------------------------------!
!   Computes the turbulent (viscosity/density) for RANS models.                !
!------------------------------------------------------------------------------!
  implicit none
!--------------------------------[Arguments]-----------------------------------!
  type(Turb_Type), target :: turb
!---------------------------------[Calling]------------------------------------!
  real :: Roughness_Coefficient
  real :: Tau_Wall_Low_Re
  real :: U_Plus_Log_Law
  real :: U_Plus_Rough_Walls
  real :: Y_Plus_Low_Re
  real :: Y_Plus_Rough_Walls
!----------------------------------[Locals]------------------------------------!
  type(Field_Type), pointer :: flow
  type(Grid_Type),  pointer :: grid
  type(Var_Type),   pointer :: u, v, w
  type(Var_Type),   pointer :: kin, eps, zeta, f22
  integer                   :: c, c1, c2, s
  real                      :: u_tan, u_tau
  real                      :: beta, pr
  real                      :: u_plus, ebf, kin_vis
  real                      :: z_o
!==============================================================================!
!   Dimensions:                                                                !
!                                                                              !
!   Production    p_kin    [m^2/s^3]   | Rate-of-strain  shear     [1/s]       !
!   Dissipation   eps % n  [m^2/s^3]   | Turb. visc.     vis_t     [kg/(m*s)]  !
!   Wall shear s. tau_wall [kg/(m*s^2)]| Dyn visc.       viscosity [kg/(m*s)]  !
!   Density       density  [kg/m^3]    | Turb. kin en.   kin % n   [m^2/s^2]   !
!   Cell volume   vol      [m^3]       | Length          lf        [m]         !
!   left hand s.  A        [kg/s]      | right hand s.   b         [kg*m^2/s^3]!
!   Wall visc.    vis_w    [kg/(m*s)]  | kinematic viscosity       [m^2/s]     !
!   Thermal cap.  capacity[m^2/(s^2*K)]| Therm. conductivity     [kg*m/(s^3*K)]!
!------------------------------------------------------------------------------!
!   p_kin = 2*vis_t / density S_ij S_ij                                        !
!   shear = sqrt(2 S_ij S_ij)                                                  !
!------------------------------------------------------------------------------!

  ! Take aliases
  flow => turb % pnt_flow
  grid => flow % pnt_grid
  call Field_Mod_Alias_Momentum   (flow, u, v, w)
  call Turb_Mod_Alias_K_Eps_Zeta_F(turb, kin, eps, zeta, f22)

  call Time_And_Length_Scale(grid, turb)

  ! Pure k-eps-zeta-f
  if(turbulence_model .eq. K_EPS_ZETA_F) then
    do c = -grid % n_bnd_cells, grid % n_cells
      turb % vis_t(c) = c_mu_d * flow % density(c) * zeta % n(c)  &
                      * kin % n(c) * turb % t_scale(c)
    end do

  ! Hybrid between k-eps-zeta-f and dynamic SGS model
  else if(turbulence_model .eq. HYBRID_LES_RANS) then
    do c = -grid % n_bnd_cells, grid % n_cells
      turb % vis_t(c) = c_mu_d * flow % density(c) * zeta % n(c)  &
                      * kin % n(c) * turb % t_scale(c)
      turb % vis_t_eff(c) = max(turb % vis_t(c),  &
                                turb % vis_t_sgs(c))
    end do
    call Comm_Mod_Exchange_Real(grid, turb % vis_t_eff)

  end if

  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    if(c2 < 0) then

      ! kinematic viscosities
      kin_vis = flow % viscosity(c1) / flow % density(c1)

      if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALL .or.  &
         Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALLFL) then

        u_tan = Field_Mod_U_Tan(flow, s)

        u_tau = c_mu25 * sqrt(kin % n(c1))
        turb % y_plus(c1) = Y_Plus_Low_Re(turb,                  &
                                          u_tau,                 &
                                          grid % wall_dist(c1),  &
                                          kin_vis)

        turb % tau_wall(c1) = Tau_Wall_Low_Re(turb,                &
                                              flow % density(c1),  &
                                              u_tau,               &
                                              u_tan,               &
                                              turb % y_plus(c1))

        ebf = 0.01 * turb % y_plus(c1) ** 4  &
            / (1.0 + 5.0 * turb % y_plus(c1))

        u_plus = U_Plus_Log_Law(turb, turb % y_plus(c1))

        if(turb % y_plus(c1) < 3.0) then
          turb % vis_w(c1) = turb % vis_t(c1) + flow % viscosity(c1)
        else
          turb % vis_w(c1) =    turb % y_plus(c1) * flow % viscosity(c1)  &
                           / (  turb % y_plus(c1) * exp(-1.0 * ebf)       &
                              + u_plus * exp(-1.0/ebf) + TINY)
        end if

        turb % y_plus(c1) = Y_Plus_Low_Re(turb,                  &
                                          u_tau,                 &
                                          grid % wall_dist(c1),  &
                                          kin_vis)

        if(turb % rough_walls) then
          z_o = Roughness_Coefficient(turb, turb % z_o_f(c1))
!         z_o = max(grid % wall_dist(c1)/(e_log * turb % y_plus(c1)),z_o) 
          turb % y_plus(c1) = Y_Plus_Rough_Walls(turb,                  &
                                                 u_tau,                 &
                                                 grid % wall_dist(c1),  &
                                                 kin_vis)
          u_plus     = U_Plus_Rough_Walls(turb, grid % wall_dist(c1))
          turb % vis_w(c1) = turb % y_plus(c1) * flow % viscosity(c1) / u_plus
        end if

        if(heat_transfer) then
          pr_t = Turb_Mod_Prandtl_Number(turb, c1)
          pr = flow % viscosity(c1) * flow % capacity / flow % conductivity
          beta = 9.24 * ((pr/pr_t)**0.75 - 1.0)     &
               * (1.0 + 0.28 * exp(-0.007*pr/pr_t))
          ebf = 0.01 * (pr * turb % y_plus(c1)**4    &
                     / ((1.0 + 5.0 * pr**3 * turb % y_plus(c1)) + TINY))
          turb % con_w(c1) =    turb % y_plus(c1)                         &
                              * flow % viscosity(c1)                      &
                              * flow % capacity                           &
                      / (  turb % y_plus(c1) * pr * exp(-1.0 * ebf)       &
                         + (u_plus + beta) * pr_t * exp(-1.0 / ebf) + TINY)
        end if
      end if  ! Grid_Mod_Bnd_Cond_Type(grid,c2).eq.WALL or WALLFL
    end if    ! c2 < 0
  end do

  call Comm_Mod_Exchange_Real(grid, turb % vis_t)
  call Comm_Mod_Exchange_Real(grid, turb % vis_w)
  if(heat_transfer) then
    call Comm_Mod_Exchange_Real(grid, turb % con_w)
  end if

  end subroutine
