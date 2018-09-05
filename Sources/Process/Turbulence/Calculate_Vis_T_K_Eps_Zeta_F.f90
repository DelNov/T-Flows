!==============================================================================!
  subroutine Calculate_Vis_T_K_Eps_Zeta_F(grid)
!------------------------------------------------------------------------------!
!   Computes the turbulent (viscosity/density) for RANS models.                !
!---------------------------------[Modules]------------------------------------!
  use Const_Mod
  use Control_Mod
  use Flow_Mod
  use Comm_Mod
  use Les_Mod
  use Rans_Mod
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!--------------------------------[Arguments]-----------------------------------!
  type(Grid_Type) :: grid
!---------------------------------[Calling]------------------------------------!
  real :: Turbulent_Prandtl_Number
!----------------------------------[Locals]------------------------------------!
  integer :: c, c1, c2, s
  real    :: beta, pr
  real    :: g_blend, y_pl, u_plus, ebf
  real    :: kin_vis
!==============================================================================!
!   Dimensions:                                                                !
!   Production    p_kin    [m^2/s^3]   | Rate-of-strain  shear     [1/s]       !
!   Dissipation   eps % n  [m^2/s^3]   | Turb. visc.     vis_t     [kg/(m*s)]  !
!   Wall shear s. tau_wall [kg/(m*s^2)]| Dyn visc.       viscosity [kg/(m*s)]  !
!   Density       density  [kg/m^3]    | Turb. kin en.   kin % n   [m^2/s^2]   !
!   Cell volume   vol      [m^3]       | Length          lf        [m]         !
!   left hand s.  A        [kg/s]      | right hand s.   b         [kg*m^2/s^3]!
!   Wall visc.    vis_wall [kg/(m*s)]  |                                       !
!   Thermal cap.  capacity[m^2/(s^2*K)]| Therm. conductivity     [kg*m/(s^3*K)]!
!------------------------------------------------------------------------------!
!   p_kin = 2*vis_t / density S_ij S_ij                                        !
!   shear = sqrt(2 S_ij S_ij)                                                  !
!------------------------------------------------------------------------------!

  call Time_And_Length_Scale(grid)

  ! c = 1, grid % n_cells
  if(turbulence_model      .eq. K_EPS_ZETA_F .and.  &
     .not. turbulence_statistics) then
    do c = -grid % n_bnd_cells, grid % n_cells
      vis_t(c) = c_mu_d * density * zeta % n(c) * kin % n(c) * t_scale(c)
    end do

  else if(turbulence_model      .eq. K_EPS_ZETA_F .and.  &
          turbulence_statistics) then
    do c = -grid % n_bnd_cells, grid % n_cells
      vis_t(c)     = c_mu_d * zeta % n(c) * kin % n(c) * t_scale(c)
      vis_t_eff(c) = max(vis_t(c), vis_t_sgs(c))
    end do
    call Comm_Mod_Exchange_Real(grid, vis_t_eff)

  end if

  ! kinematic viscosities
  kin_vis = viscosity / density

  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    if(c2 < 0 .and. Grid_Mod_Bnd_Cond_Type(grid,c2) .ne. BUFFER) then
      if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALL .or.  &
         Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALLFL) then

        u_tau(c1)  = c_mu**0.25 * kin % n(c1)**0.5

        if(rough_walls) then
          y_plus(c1) = (grid % wall_dist(c1)+Zo)*u_tau(c1)/kin_vis
        else if(.not. rough_walls) then
          y_plus(c1) = grid % wall_dist(c1)*u_tau(c1)/kin_vis
        end if

        g_blend  = 0.01*y_plus(c1)**4 / (1.0 + 5.0*y_plus(c1))

        y_pl   = max(y_plus(c1), 0.13)
        u_plus = log(y_pl*e_log)/kappa

        if(y_pl < 3.0) then
          vis_wall(c1) = vis_t(c1) + viscosity
        else
          vis_wall(c1) = y_plus(c1)*viscosity/ &
            (y_pl*exp(-1.0*g_blend) + u_plus*exp(-1.0/g_blend) + TINY)
        end if

        if(rough_walls) then
          u_plus = log((grid % wall_dist(c1)+Zo)/Zo)/(kappa + TINY) + TINY
          vis_wall(c1) = min(y_pl*viscosity*kappa/ &
            log((grid % wall_dist(c1)+Zo)/Zo), MEGA*kin_vis)
        end if

        if(heat_transfer) then
          pr_t = Turbulent_Prandtl_Number(grid, c1)
          pr = viscosity * capacity / conductivity
          beta = 9.24 * ((pr/pr_t)**0.75 - 1.0) * &
            (1.0 + 0.28 * exp(-0.007*pr/pr_t))
          ebf = 0.01 * (pr*y_pl)**4 / &
            ((1.0 + 5.0 * pr**3 * y_pl) + TINY)
          con_wall(c1) = y_pl*viscosity*capacity/(y_pl*pr* &
            exp(-1.0 * ebf) + (u_plus + beta)*pr_t*exp(-1.0/ebf) + TINY)
        end if
      end if  ! Grid_Mod_Bnd_Cond_Type(grid,c2).eq.WALL or WALLFL
    end if    ! c2 < 0
  end do

  call Comm_Mod_Exchange_Real(grid, vis_t)
  call Comm_Mod_Exchange_Real(grid, vis_wall)

  end subroutine
