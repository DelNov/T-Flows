!==============================================================================!
  subroutine Vis_T_K_Eps_Zeta_F(Turb)
!------------------------------------------------------------------------------!
!   Computes the turbulent (viscosity/density) for RANS models.                !
!------------------------------------------------------------------------------!
  implicit none
!--------------------------------[Arguments]-----------------------------------!
  class(Turb_Type), target :: Turb
!----------------------------------[Locals]------------------------------------!
  type(Field_Type), pointer :: Flow
  type(Grid_Type),  pointer :: Grid
  type(Var_Type),   pointer :: u, v, w
  type(Var_Type),   pointer :: kin, eps, zeta, f22
  integer                   :: c, c1, c2, s, reg
  real                      :: u_tan, u_tau
  real                      :: beta, pr, sc
  real                      :: u_plus, ebf, kin_vis
  real                      :: z_o
!------------------------------------------------------------------------------!
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
!==============================================================================!

  ! Take aliases
  Flow => Turb % pnt_flow
  Grid => Flow % pnt_grid
  call Flow % Alias_Momentum    (u, v, w)
  call Turb % Alias_K_Eps_Zeta_F(kin, eps, zeta, f22)

  call Turb % Time_And_Length_Scale(Grid)

  ! Pure k-eps-zeta-f
  if(Turb % model .eq. K_EPS_ZETA_F) then
    do reg = Boundary_And_Inside_Regions()
      do c = Cells_In_Region(reg)
        Turb % vis_t(c) = Turb % c_mu_d * Flow % density(c) * zeta % n(c)  &
                        * kin % n(c) * Turb % t_scale(c)
      end do
    end do
    call Grid % Exchange_Cells_Real(Turb % vis_t)

  ! Hybrid between k-eps-zeta-f and dynamic SGS model
  else if(Turb % model .eq. HYBRID_LES_RANS) then
    do reg = Boundary_And_Inside_Regions()
      do c = Cells_In_Region(reg)
        Turb % vis_t(c) = Turb % c_mu_d * Flow % density(c) * zeta % n(c)  &
                        * kin % n(c) * Turb % t_scale(c)
        Turb % vis_t_eff(c) = max(Turb % vis_t(c),  &
                                  Turb % vis_t_sgs(c))
      end do
    end do
    call Grid % Exchange_Cells_Real(Turb % vis_t)
    call Grid % Exchange_Cells_Real(Turb % vis_t_eff)

  end if

  do reg = Boundary_Regions()
    if(Grid % region % type(reg) .eq. WALL .or.  &
       Grid % region % type(reg) .eq. WALLFL) then
      do s = Faces_In_Region(reg)
        c1 = Grid % faces_c(1,s)
        c2 = Grid % faces_c(2,s)

        Assert(c2 < 0)

        kin_vis =  Flow % viscosity(c1) / Flow % density(c1)

        ! Set up roughness coefficient
        z_o = Turb % Roughness_Coeff(c1, c2)

        ! Compute tangential velocity component
        u_tan = Flow % U_Tan(s)

        u_tau = Turb % c_mu25 * sqrt(kin % n(c1))

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
          pr_t = Turb % Prandtl_Turb(c1)
          pr   = Flow % Prandtl_Numb(c1)          ! laminar Prandtl number
          beta = Turb % Beta_Scalar(pr, pr_t)
          ! According to Toparlar et al. 2019 paper
          ! "CFD simulation of the near-neutral atmospheric boundary layer:
          ! New temperature inlet profile consistent with wall functions"

          if(z_o .gt. TINY) then
            beta = 0.0
          end if

          ebf = Turb % Ebf_Scalar(c1, pr)

          Turb % con_w(c1) =    Turb % y_plus(c1)                           &
                              * Flow % viscosity(c1)                        &
                              * Flow % capacity(c1)                         &
                      / ( Turb % y_plus(c1) * pr * exp(-1.0 * ebf)          &
                        + (u_plus + beta) * pr_t * exp(-1.0 / ebf) + TINY )
        end if

        if(Flow % n_scalars .gt. 0) then
          sc   = Flow % Schmidt_Numb(c1)          ! laminar Schmidt number
          beta = Turb % Beta_Scalar(sc, sc_t)
          ! According to Toparlar et al. 2019 paper
          ! "CFD simulation of the near-neutral atmospheric boundary layer:
          ! New temperature inlet profile consistent with wall functions"

          if(z_o .gt. TINY) then
            beta = 0.0
          end if

          ebf = Turb % Ebf_Scalar(c1, pr)

          Turb % diff_w(c1) =  Turb % y_plus(c1)                  &
              * (Flow % viscosity(c1)/Flow % density(c1))         &
              / (Turb % y_plus(c1) * sc * exp(-1.0 * ebf)         &
               + (u_plus + beta) * sc_t * exp(-1.0 / ebf) + TINY)
        end if

      end do    ! faces in regions
    end if      ! region is WALL or WALLFL
  end do        ! through regions

  call Grid % Exchange_Cells_Real(Turb % vis_w)
  if(Flow % heat_transfer) then
    call Grid % Exchange_Cells_Real(Turb % con_w)
  end if
  if(Flow % n_scalars > 0) then
    call Grid % Exchange_Cells_Real(Turb % diff_w)
  end if

  end subroutine
