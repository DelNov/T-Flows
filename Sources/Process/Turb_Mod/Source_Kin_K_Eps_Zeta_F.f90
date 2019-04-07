!==============================================================================!
  subroutine Source_Kin_K_Eps_Zeta_F(flow, sol)
!------------------------------------------------------------------------------!
!   Computes the source terms in kin transport equation.                       !
!------------------------------------------------------------------------------!
!   In kinetic energy equation there are two source terms:                     !
!                                                                              !
!     /
!    |                                                                         !
!    | (density (p_kin - eps)) dV                                              !
!    |                                                                         !
!   /                                                                          !
!                                                                              !
!------------------------------------------------------------------------------!
!---------------------------------[Modules]------------------------------------!
  use Const_Mod
  use Field_Mod
  use Turbulence_Mod
  use Grid_Mod,       only: Grid_Type
  use Solver_Mod,     only: Solver_Type
  use Matrix_Mod,     only: Matrix_Type
!------------------------------------------------------------------------------!
  implicit none
!--------------------------------[Arguments]-----------------------------------!
  type(Field_Type),  target :: flow
  type(Solver_Type), target :: sol
!---------------------------------[Calling]------------------------------------!
  real :: Y_Plus_Low_Re
  real :: Y_Plus_Rough_Walls
  real :: Roughness_Coefficient
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),   pointer :: grid
  type(Var_Type),    pointer :: u, v, w
  type(Matrix_Type), pointer :: a
  real,              pointer :: b(:)
  integer                    :: c, c1, c2, s
  real                       :: u_tan, u_nor_sq, u_nor, u_tot_sq
  real                       :: lf, ebf, p_kin_int, p_kin_wf
  real                       :: alpha1, l_rans, l_sgs, kin_vis
!==============================================================================!
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
!------------------------------------------------------------------------------!

  ! Take aliases
  grid => flow % pnt_grid
  u    => flow % u
  v    => flow % v
  w    => flow % w
  a    => sol % a
  b    => sol % b % val

  ! Production source:
  do c = 1, grid % n_cells
    p_kin(c) = vis_t(c) * shear(c)**2
    b(c)     = b(c) + p_kin(c) * grid % vol(c)
  end do

  if(turbulence_model .eq. HYBRID_LES_RANS) then
    do c = 1, grid % n_cells
      lf = grid % vol(c)**ONE_THIRD
      l_sgs  = 0.8*lf
      l_rans = 0.41*grid % wall_dist(c)
      alpha1 = max(1.0,l_rans/l_sgs)

      if(alpha1 < 1.05) then
        a % val(a % dia(c)) = a % val(a % dia(c))   &
                            + density * eps % n(c)  &
                            / (kin % n(c) + TINY) * grid % vol(c)
      else
        a % val(a % dia(c)) = a % val(a % dia(c))   &
          + density                                 &
          * min(alpha1**1.45 * eps % n(c), kin % n(c)**1.5 / (lf*0.01))  &
          / (kin % n(c) + TINY) * grid % vol(c)
      end if
    end do
  else  ! turbuence model will be K_EPS_ZETA_F
    do c = 1, grid % n_cells
      a % val(a % dia(c)) = a % val(a % dia(c))   &
                          + density * eps % n(c)  &
                          / (kin % n(c) + TINY) * grid % vol(c)

      if(buoyancy) then
        buoy_beta(c) = 1.0
        g_buoy(c) = -buoy_beta(c) * (grav_x * ut % n(c) +  &
                                     grav_y * vt % n(c) +  &
                                     grav_z * wt % n(c)) * density
        b(c) = b(c) + max(0.0, g_buoy(c) * grid % vol(c))
        a % val(a % dia(c)) = a % val(a % dia(c))  &
                  + max(0.0,-g_buoy(c) * grid % vol(c) / (kin % n(c) + TINY))
      end if
    end do
  end if

  ! Kinematic viscosities
  kin_vis = viscosity / density

  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    if(c2 < 0) then
      if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALL .or. &
         Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALLFL) then

        ! Compute tangential velocity component
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

        if(rough_walls) then
          z_o = Roughness_Coefficient(grid, z_o_f(c1), c1)       
          u_tau(c1)  = c_mu25 * sqrt(kin % n(c1))
          y_plus(c1) = Y_Plus_Rough_Walls(u_tau(c1), &
                       grid % wall_dist(c1), kin_vis) 

          tau_wall(c1) = density*kappa*u_tau(c1)*u_tan  &
                       / log(((grid % wall_dist(c1)+z_o) / z_o))

          p_kin(c1) = tau_wall(c1) * c_mu25 * sqrt(kin % n(c1)) &
                      / (kappa*(grid % wall_dist(c1)+z_o))
          b(c1)     = b(c1) + (p_kin(c1)  &
                    - vis_t(c1) * shear(c1)**2) * grid % vol(c1)
        else
          u_tau(c1) = c_mu25 * sqrt(kin % n(c1))
          y_plus(c1) = Y_Plus_Low_Re(u_tau(c1), grid % wall_dist(c1), kin_vis)

          tau_wall(c1) = density*kappa*u_tau(c1)*u_tan  &
                       / log(e_log*max(y_plus(c1),1.05))

          ebf = max(0.01 * y_plus(c1)**4 / (1.0 + 5.0*y_plus(c1)),tiny)

          p_kin_wf  = tau_wall(c1) * 0.07**0.25 * sqrt(kin % n(c1))  &
                    / (grid % wall_dist(c1) * kappa)

          p_kin_int = vis_t(c1) * shear(c1)**2

          p_kin(c1) = p_kin_int * exp(-1.0 * ebf) + p_kin_wf  &
                    * exp(-1.0 / ebf)
          b(c1)     = b(c1) + (p_kin(c1) - p_kin_int) * grid % vol(c1)
        end if! rough_walls
      end if  ! Grid_Mod_Bnd_Cond_Type(grid,c2).eq.WALL or WALLFL
    end if    ! c2 < 0
  end do

  end subroutine
