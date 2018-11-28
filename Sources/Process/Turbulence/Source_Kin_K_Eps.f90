!==============================================================================!
  subroutine Source_Kin_K_Eps(grid, sol)
!------------------------------------------------------------------------------!
!   Computes the source terms in kin transport equation for k-epsilon model    !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
  use Flow_Mod
  use Comm_Mod
  use Les_Mod
  use Rans_Mod
  use Grad_Mod
  use Grid_Mod,   only: Grid_Type
  use Solver_Mod, only: Solver_Type
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type)           :: grid
  type(Solver_Type), target :: sol
!---------------------------------[Calling]------------------------------------!
  real :: Y_Plus_Low_Re
  real :: Roughness_Coefficient
!-----------------------------------[Locals]-----------------------------------!
  type(Matrix_Type), pointer :: a
  real,              pointer :: b(:)
  integer                    :: c, c1, c2, s
  real                       :: u_tot2, u_nor, u_nor2, u_tan
  real                       :: kin_vis  ! [m^2/s]
  real                       :: ebf, p_kin_int, p_kin_wf, u_tau_new
!==============================================================================!
!   Dimensions:                                                                !
!                                                                              !
!   production    p_kin    [m^2/s^3]   | rate-of-strain  shear     [1/s]       !
!   dissipation   eps % n  [m^2/s^3]   | turb. visc.     vis_t     [kg/(m*s)]  !
!   wall shear s. tau_wall [kg/(m*s^2)]| dyn visc.       viscosity [kg/(m*s)]  !
!   density       density  [kg/m^3]    | turb. kin en.   kin % n   [m^2/s^2]   !
!   cell volume   vol      [m^3]       | length          lf        [m]         !
!   left hand s.  a        [kg/s]      | right hand s.   b         [kg*m^2/s^3]!
!------------------------------------------------------------------------------!
!   p_kin = 2*vis_t / density S_ij S_ij                                        !
!   shear = sqrt(2 S_ij S_ij)                                                  !
!------------------------------------------------------------------------------!

  ! Take aliases
  a => sol % a
  b => sol % b

  !-----------------------------------------!
  !   Compute the sources in the interior   !
  !-----------------------------------------!
  ! Production source:
  do c = 1, grid % n_cells
    p_kin(c) = vis_t(c) * shear(c)**2
    b(c)     = b(c) + p_kin(c) * grid % vol(c)

    a % val(a % dia(c)) = a % val(a % dia(c)) + &
         density * eps % n(c)/(kin % n(c) + TINY) * grid % vol(c)

    if (buoyancy) then
      buoy_beta(c) = 1.0
      g_buoy(c) = -buoy_beta(c) * (grav_x * ut % n(c) +  &
                                   grav_y * vt % n(c) +  &
                                   grav_z * wt % n(c)) * density
      b(c) = b(c) + max(0.0, g_buoy(c) * grid % vol(c))
      a % val(a % dia(c)) = a % val(a % dia(c))  &
      + max(0.0,-g_buoy(c) * grid % vol(c) / (kin % n(c) + TINY))
    end if
  end do

  kin_vis = viscosity / density

  !-----------------------------------------------!
  !  Compute the sources in the near wall cells   !
  !-----------------------------------------------!
  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    if(c2 < 0) then
      if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALL .or.  &
         Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALLFL) then

        ! Compute tangential velocity component
        u_tot2 = u % n(c1) **2  &
               + v % n(c1) **2  &
               + w % n(c1) **2
        u_nor = ( u % n(c1) * grid % sx(s)     &
                + v % n(c1) * grid % sy(s)     &
                + w % n(c1) * grid % sz(s) )   &
                     / sqrt(  grid % sx(s)**2  &
                            + grid % sy(s)**2  &
                            + grid % sz(s)**2 )
        u_nor2 = u_nor**2

        if(u_tot2  >  u_nor2) then
          u_tan = sqrt(u_tot2 - u_nor2)
        else
          u_tan = TINY
        end if

        if(rough_walls) then
          z_o = Roughness_Coefficient(grid, c1)      
          u_tau(c1) = c_mu25 * sqrt(kin % n(c1))
          y_plus(c1) = u_tau(c1) * (grid % wall_dist(c1) + z_o) &
                     / kin_vis

          tau_wall(c1) = density*kappa*u_tau(c1)*u_tan  &
                       / log(max(((grid % wall_dist(c1)+z_o)/z_o), 1.05))

          p_kin(c1) = density * tau_wall(c1) * c_mu25 * sqrt(kin % n(c1))   &
                    / (kappa*(grid % wall_dist(c1) + z_o))
          b(c1)     = b(c1) + (p_kin(c1) - p_kin(c1))   &
                    * grid % vol(c1)
        else
          u_tau(c1) = c_mu25 * sqrt(kin % n(c1))
          y_plus(c1) = Y_Plus_Low_Re(u_tau(c1), grid % wall_dist(c1), kin_vis)

          tau_wall(c1) = density*kappa*u_tau(c1)*u_tan   &
                       / log(e_log * max(y_plus(c1), 1.05))

          ebf = 0.01 * y_plus(c1)**4 / (1.0 + 5.0*y_plus(c1))

          p_kin_wf  = tau_wall(c1) * 0.07**0.25 * sqrt(kin % n(c1))  &
                    / (grid % wall_dist(c1) * kappa)

          p_kin_int = vis_t(c1) * shear(c1)**2

          p_kin(c1) = p_kin_int * exp(-1.0 * ebf) + p_kin_wf  &
                    * exp(-1.0 / ebf)

          b(c1) = b(c1) + (p_kin(c1) - p_kin_int) * grid % vol(c1)
        end if  ! rough_walls
      end if    ! Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALL or WALLFL
    end if      ! c2 < 0
  end do

  call Comm_Mod_Exchange_Real(grid, kin % n)

  end subroutine
