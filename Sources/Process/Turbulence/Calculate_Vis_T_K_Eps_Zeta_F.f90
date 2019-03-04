!==============================================================================!
  subroutine Calculate_Vis_T_K_Eps_Zeta_F(flow)
!------------------------------------------------------------------------------!
!   Computes the turbulent (viscosity/density) for RANS models.                !
!---------------------------------[Modules]------------------------------------!
  use Const_Mod
  use Control_Mod
  use Field_Mod
  use Comm_Mod
  use Les_Mod
  use Rans_Mod
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!--------------------------------[Arguments]-----------------------------------!
  type(Field_Type), target :: flow
!---------------------------------[Calling]------------------------------------!
  real :: Turbulent_Prandtl_Number
  real :: U_Plus_Log_Law
  real :: U_Plus_Rough_Walls
  real :: Y_Plus_Low_Re
  real :: Y_Plus_Rough_Walls
  real :: Roughness_Coefficient
!----------------------------------[Locals]------------------------------------!
  type(Grid_Type), pointer :: grid
  type(Var_Type),  pointer :: u, v, w
  integer                  :: c, c1, c2, s
  real                     :: u_tan, u_nor_sq, u_nor, u_tot_sq
  real                     :: beta, pr
  real                     :: u_plus, ebf, kin_vis
!==============================================================================!
!   Dimensions:                                                                !
!                                                                              !
!   Production    p_kin    [m^2/s^3]   | Rate-of-strain  shear     [1/s]       !
!   Dissipation   eps % n  [m^2/s^3]   | Turb. visc.     vis_t     [kg/(m*s)]  !
!   Wall shear s. tau_wall [kg/(m*s^2)]| Dyn visc.       viscosity [kg/(m*s)]  !
!   Density       density  [kg/m^3]    | Turb. kin en.   kin % n   [m^2/s^2]   !
!   Cell volume   vol      [m^3]       | Length          lf        [m]         !
!   left hand s.  A        [kg/s]      | right hand s.   b         [kg*m^2/s^3]!
!   Wall visc.    vis_wall [kg/(m*s)]  | kinematic viscosity       [m^2/s]     !
!   Thermal cap.  capacity[m^2/(s^2*K)]| Therm. conductivity     [kg*m/(s^3*K)]!
!------------------------------------------------------------------------------!
!   p_kin = 2*vis_t / density S_ij S_ij                                        !
!   shear = sqrt(2 S_ij S_ij)                                                  !
!------------------------------------------------------------------------------!

  ! Take aliases
  grid => flow % pnt_grid
  u    => flow % u
  v    => flow % v
  w    => flow % w

  call Time_And_Length_Scale(grid)

  ! Pure k-eps-zeta-f
  if(turbulence_model .eq. K_EPS_ZETA_F) then
    do c = -grid % n_bnd_cells, grid % n_cells
      vis_t(c) = c_mu_d * density * zeta % n(c) * kin % n(c) * t_scale(c)
    end do

  ! Hybrid between k-eps-zeta-f and dynamic SGS model
  else if(turbulence_model .eq. HYBRID_LES_RANS) then
    do c = -grid % n_bnd_cells, grid % n_cells
      vis_t(c)     = c_mu_d * density * zeta % n(c)  &
                   * kin % n(c) * t_scale(c)
      vis_t_eff(c) = max(vis_t(c), vis_t_sgs(c))
    end do
    call Comm_Mod_Exchange_Real(grid, vis_t_eff)

  end if

  ! kinematic viscosities
  kin_vis = viscosity / density

  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    if(c2 < 0) then
      if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALL .or.  &
         Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALLFL) then

        u_tot_sq = u % n(c1) * u % n(c1) &
                 + v % n(c1) * v % n(c1) &
                 + w % n(c1) * w % n(c1)
        u_nor  = ( u % n(c1) * grid % sx(s)     &
                 + v % n(c1) * grid % sy(s)     &
                 + w % n(c1) * grid % sz(s) )   &
                 / sqrt(  grid % sx(s)*grid % sx(s)  &
                        + grid % sy(s)*grid % sy(s)  &
                        + grid % sz(s)*grid % sz(s))
        u_nor_sq = u_nor**2

        if( u_tot_sq  > u_nor_sq) then
          u_tan = sqrt(u_tot_sq - u_nor_sq)
        else
          u_tan = TINY
        end if

        u_tau(c1) = c_mu25 * sqrt(kin % n(c1))
        y_plus(c1) = Y_Plus_Low_Re(u_tau(c1), grid % wall_dist(c1), kin_vis)

        tau_wall(c1) = density*kappa*u_tau(c1)*u_tan    &
                     / log(e_log*max(y_plus(c1),1.05))

        ebf = 0.01 * y_plus(c1)**4 / (1.0 + 5.0*y_plus(c1))

        u_plus = U_Plus_Log_Law(y_plus(c1))

        if(y_plus(c1) < 3.0) then
          vis_wall(c1) = vis_t(c1) + viscosity
        else
          vis_wall(c1) = y_plus(c1) * viscosity         &
                       / (  y_plus(c1) * exp(-1.0*ebf)  &
                          + u_plus     * exp(-1.0/ebf) + TINY)
        end if

        y_plus(c1) = Y_Plus_Low_Re(u_tau(c1), grid % wall_dist(c1), kin_vis)

        if(rough_walls) then
          z_o = Roughness_Coefficient(grid, z_o_f(c1), c1)      
          y_plus(c1) = Y_Plus_Rough_Walls(u_tau(c1),             &
                                          grid % wall_dist(c1),  &
                                          kin_vis)
          u_plus     = U_Plus_Rough_Walls(grid % wall_dist(c1))
          vis_wall(c1) = y_plus(c1) * viscosity * kappa  &
                       / log((grid % wall_dist(c1)+z_o)/z_o)  ! is this U+?
        end if  

        if(heat_transfer) then
          pr_t = Turbulent_Prandtl_Number(grid, c1)
          pr = viscosity * capacity / conductivity
          beta = 9.24 * ((pr/pr_t)**0.75 - 1.0) * &
            (1.0 + 0.28 * exp(-0.007*pr/pr_t))
          ebf = 0.01 * (pr*y_plus(c1)**4 / &
            ((1.0 + 5.0 * pr**3 * y_plus(c1)) + TINY))
          con_wall(c1) =    y_plus(c1) * viscosity * capacity          &
                       / (  y_plus(c1) * pr        * exp(-1.0 * ebf)   &
                          +(u_plus + beta) * pr_t  * exp(-1.0/ebf) + TINY)
        end if
      end if  ! Grid_Mod_Bnd_Cond_Type(grid,c2).eq.WALL or WALLFL
    end if    ! c2 < 0
  end do

  call Comm_Mod_Exchange_Real(grid, vis_t)
  call Comm_Mod_Exchange_Real(grid, vis_wall)
  if(heat_transfer) then
    call Comm_Mod_Exchange_Real(grid, con_wall)
  end if

  end subroutine
