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
  type(Field_Type), pointer :: Flow
  type(Grid_Type),  pointer :: Grid
  type(Var_Type),   pointer :: u, v, w
  type(Var_Type),   pointer :: kin, eps, zeta, f22
  integer                   :: c, c1, c2, s
  real                      :: u_tan, u_tau
  real                      :: beta, pr, sc
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
  Flow => turb % pnt_flow
  Grid => Flow % pnt_grid
  call Flow % Alias_Momentum(u, v, w)
  call Turb_Mod_Alias_K_Eps_Zeta_F(turb, kin, eps, zeta, f22)

  call Time_And_Length_Scale(Grid, turb)

  ! Pure k-eps-zeta-f
  if(turb % model .eq. K_EPS_ZETA_F) then
    do c = -Grid % n_bnd_cells, Grid % n_cells
      turb % vis_t(c) = c_mu_d * Flow % density(c) * zeta % n(c)  &
                      * kin % n(c) * turb % t_scale(c)
    end do

  ! Hybrid between k-eps-zeta-f and dynamic SGS model
  else if(turb % model .eq. HYBRID_LES_RANS) then
    do c = -Grid % n_bnd_cells, Grid % n_cells
      turb % vis_t(c) = c_mu_d * Flow % density(c) * zeta % n(c)  &
                      * kin % n(c) * turb % t_scale(c)
      turb % vis_t_eff(c) = max(turb % vis_t(c),  &
                                turb % vis_t_sgs(c))
    end do
    call Grid % Exchange_Cells_Real(turb % vis_t_eff)

  end if

  do s = 1, Grid % n_faces
    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)

    if(c2 < 0) then

      ! kinematic viscosities
      kin_vis = Flow % viscosity(c1) / Flow % density(c1)

      if(Grid % Bnd_Cond_Type(c2) .eq. WALL .or.  &
         Grid % Bnd_Cond_Type(c2) .eq. WALLFL) then

        u_tan = Flow % U_Tan(s)

        u_tau = c_mu25 * sqrt(kin % n(c1))
        turb % y_plus(c1) = Y_Plus_Low_Re(turb,                  &
                                          u_tau,                 &
                                          Grid % wall_dist(c1),  &
                                          kin_vis)

        turb % tau_wall(c1) = Tau_Wall_Low_Re(turb,                &
                                              Flow % density(c1),  &
                                              u_tau,               &
                                              u_tan,               &
                                              turb % y_plus(c1))

        ebf = Turb_Mod_Ebf_Momentum(turb, c1)

        u_plus = U_Plus_Log_Law(turb, turb % y_plus(c1))

        if(turb % y_plus(c1) < 3.0) then
          turb % vis_w(c1) = turb % vis_t(c1) + Flow % viscosity(c1)
        else

          if(turb % rough_walls) then
            z_o = Roughness_Coefficient(turb, turb % z_o_f(c1))
            z_o = max(Grid % wall_dist(c1)   &
                / (e_log * max(turb % y_plus(c1),1.0)), z_o)

            turb % y_plus(c1) = Y_Plus_Rough_Walls(turb,                  &
                                                   u_tau,                 &
                                                   Grid % wall_dist(c1),  &
                                                   kin_vis)
            u_plus     = U_Plus_Rough_Walls(turb, Grid % wall_dist(c1))
          end if

          if(turb % y_plus(c1) < 11.3 ) then
            ebf = 0.00001 
          end if
         
          turb % vis_w(c1) =    turb % y_plus(c1) * Flow % viscosity(c1)  &
                           / (  turb % y_plus(c1) * exp(-1.0 * ebf)       &
                              + u_plus * exp(-1.0/ebf) + TINY)

        end if

        if(Flow % heat_transfer) then
          pr_t = Turb_Mod_Prandtl_Number(turb, c1)
          pr   = Flow % Prandtl_Number(c1)          ! laminar Prandtl number
          beta = 9.24 * ((pr/pr_t)**0.75 - 1.0)     &
               * (1.0 + 0.28 * exp(-0.007*pr/pr_t))

          ! According to Toparlar et al. 2019 paper
          ! "CFD simulation of the near-neutral atmospheric boundary layer: New
          ! temperature inlet profile consistent with wall functions"

          if(turb % rough_walls) beta = 0.0         

          ebf = Turb_Mod_Ebf_Scalar(turb, c1, pr)
          
          turb % con_w(c1) =    turb % y_plus(c1)                           &
                              * Flow % viscosity(c1)                        &
                              * Flow % capacity(c1)                         &
                      / ( turb % y_plus(c1) * pr * exp(-1.0 * ebf)          &
                         + (u_plus + beta) * pr_t * exp(-1.0 / ebf) + TINY )
        end if

        if(Flow % n_scalars > 0) then
          sc   = Flow % Schmidt_Number(c1)          ! laminar Schmidt number
          beta = 9.24 * ((sc/sc_t)**0.75 - 1.0)                   &
               * (1.0 + 0.28 * exp(-0.007*sc/sc_t))

          ! According to Toparlar et al. 2019 paper
          ! "CFD simulation of the near-neutral atmospheric boundary layer: New
          ! temperature inlet profile consistent with wall functions"

          if(turb % rough_walls) beta = 0.0         

          ebf = 0.01 * (sc * turb % y_plus(c1)**4                 &
              / ((1.0 + 5.0 * sc**3 * turb % y_plus(c1)) + TINY))
          turb % diff_w(c1) =  turb % y_plus(c1)                  &
              * (Flow % viscosity(c1)/Flow % density(c1))         &
              / (turb % y_plus(c1) * sc * exp(-1.0 * ebf)         &
              + (u_plus + beta) * sc_t * exp(-1.0 / ebf) + TINY)
        end if

      end if  ! Grid % Bnd_Cond_Type(c2).eq.WALL or WALLFL
    end if    ! c2 < 0
  end do

  call Grid % Exchange_Cells_Real(turb % vis_t)
  call Grid % Exchange_Cells_Real(turb % vis_w)
  if(Flow % heat_transfer) then
    call Grid % Exchange_Cells_Real(turb % con_w)
  end if
  if(Flow % n_scalars > 0) then
    call Grid % Exchange_Cells_Real(turb % diff_w)
  end if

  end subroutine
