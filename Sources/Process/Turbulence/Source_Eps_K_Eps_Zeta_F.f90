!==============================================================================!
  subroutine Source_Eps_K_Eps_Zeta_F(flow, sol)
!------------------------------------------------------------------------------!
!   Calculates source terms in equation of dissipation of turbulent energy     !
!   and imposes boundary condition                                             !
!------------------------------------------------------------------------------!
!---------------------------------[Modules]------------------------------------!
  use Const_Mod
  use Grid_Mod,   only: Grid_Type
  use Field_Mod,  !only: Field_Type, viscosity, density, buoyancy
  use Solver_Mod, only: Solver_Type
  use Matrix_Mod, only: Matrix_Type
  use Les_Mod
  use Rans_Mod
!------------------------------------------------------------------------------!
  implicit none
!--------------------------------[Arguments]-----------------------------------!
  type(Field_Type),  target :: flow
  type(Solver_Type), target :: sol
!---------------------------------[Calling]------------------------------------!
  real :: Y_Plus_Low_Re
  real :: Roughness_Coefficient
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),   pointer :: grid
  type(Var_Type),    pointer :: u, v, w, t
  type(Matrix_Type), pointer :: a
  real,              pointer :: b(:)
  integer                    :: c, s, c1, c2, j
  real                       :: u_tan, u_nor_sq, u_nor, u_tot_sq
  real                       :: e_sor, c_11e
  real                       :: eps_wf, eps_int, p_kin_int, p_kin_wf, ebf
  real                       :: fa, u_tau_new, kin_vis, g_buoy_wall
!==============================================================================!
!   In dissipation of turbulent kinetic energy equation exist two              !
!   source terms which have form:                                              !
!                                                                              !
!    int( density ((Cv_e1 * p_kin - Cv_11 eps) / t_scale) * dV                 !
!                                                                              !
!   First, positive , source term is solved and added to source  coefficient   !
!   b(c) on right hand side.  Second, negative, source term is added to main   !
!   diagonal left hand side coefficient matrix in order to increase stability  !
!   of solver.  It is nessesary to calculate coefficient Cv_11 using kin,      !
!   Cv_e2, vi2 and coefficient A1                                              !
!------------------------------------------------------------------------------!
!   Dimensions:                                                                !
!                                                                              !
!   production    p_kin    [m^2/s^3]   | rate-of-strain  shear     [1/s]       !
!   dissipation   eps % n  [m^2/s^3]   | turb. visc.     vis_t     [kg/(m*s)]  !
!   wall shear s. tau_wall [kg/(m*s^2)]| dyn visc.       viscosity [kg/(m*s)]  !
!   density       density  [kg/m^3]    | turb. kin en.   kin % n   [m^2/s^2]   !
!   cell volume   vol      [m^3]       | length          lf        [m]         !
!   left hand s.  a        [kg/s]      | right hand s.   b         [kg*m^2/s^4]!
!------------------------------------------------------------------------------!
!   p_kin = 2*vis_t / density S_ij S_ij                                        !
!   shear = sqrt(2 S_ij S_ij)                                                  !
!------------------------------------------------------------------------------!

  ! Take aliases
  grid => flow % pnt_grid
  u    => flow % u
  v    => flow % v
  w    => flow % w
  t    => flow % t
  a    => sol % a
  b    => sol % b % val

  call Time_And_Length_Scale(grid)

  do c = 1, grid % n_cells 
    e_sor = grid % vol(c)/(t_scale(c)+TINY)
    c_11e = c_1e*(1.0 + alpha * ( 1.0/(zeta % n(c)+TINY) ))
    b(c) = b(c) + c_11e * e_sor * p_kin(c)

    ! Fill in a diagonal of coefficient matrix
    a % val(a % dia(c)) =  a % val(a % dia(c)) + c_2e * e_sor * density

    ! Add buoyancy (linearly split) to eps equation as required in the t2 model
    if(buoyancy) then
      b(c) = b(c) + max(0.0, c_11e * e_sor * g_buoy(c))
      a % val(a % dia(c)) = a % val(a % dia(c))  &
              + max(0.0,-c_11e * e_sor * g_buoy(c) / (eps % n(c) + TINY))
    end if
  end do

  !-------------------------------------------------------!
  !   Following block shows density dependent behaviour   !
  !-------------------------------------------------------!

  kin_vis = viscosity / density

  ! Imposing a boundary condition on wall for eps
  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    if(c2 < 0) then
      if( Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALL .or.  &
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
          eps % n(c1) = c_mu75 * kin % n(c1)**1.5 / & 
                      ((grid % wall_dist(c1) + z_o) * kappa)

          ! Adjusting coefficient to fix eps value in near wall calls
          do j = a % row(c1), a % row(c1 + 1) - 1 
            a % val(j) = 0.0 
          end do

          b(c1) = eps % n(c1) 
          a % val(a % dia(c1)) = 1.0 
        else
          u_tau(c1) = c_mu25 * sqrt(kin % n(c1))
          y_plus(c1) = Y_Plus_Low_Re(u_tau(c1), grid % wall_dist(c1), kin_vis)

          tau_wall(c1) = density * kappa * u_tau(c1) * u_tan  &
                       / log(e_log * max(y_plus(c1), 1.05))

          u_tau_new = sqrt(tau_wall(c1)/density)

          y_plus(c1) = Y_Plus_Low_Re(u_tau_new, grid % wall_dist(c1), kin_vis)

          eps_int = 2.0* kin_vis * kin % n(c1)  &
                  / grid % wall_dist(c1)**2

          eps_wf  = c_mu75 * kin % n(c1)**1.5            &
                  / (grid % wall_dist(c1) * kappa)

          ebf = max(0.01 * y_plus(c1)**4 / (1.0 + 5.0*y_plus(c1)), TINY)

          p_kin_wf  = tau_wall(c1) * c_mu25 * sqrt(kin % n(c1))  &
                      / (grid % wall_dist(c1) * kappa)
 
          p_kin_int = vis_t(c1) * shear(c1)**2
 
          p_kin(c1) = p_kin_int * exp(-1.0 * ebf) + p_kin_wf  &
                      * exp(-1.0 / ebf)

          fa = min(p_kin_wf * exp(-1.0 / ebf)/p_kin(c1) , 1.0) 

          if(buoyancy) then
            g_buoy_wall = beta_tec*abs(grav_x + grav_y + grav_z)  &
                         * sqrt(abs(t % q(c2))                    &
                         * c_mu_theta5                            &
                         * sqrt(abs(t2 % n(c1) * kin % n(c1))))

            fa = min((p_kin_wf + g_buoy_wall) * exp(-1.0 / ebf)     &
                 /(p_kin(c1) + g_buoy(c1)) , 1.0) 

          end if        

          if(y_plus(c1) > 3) then

            eps % n(c1) = sqrt(1.0 - fa) * eps_int + sqrt(fa) * eps_wf

            ! Adjusting coefficient to fix eps value in near wall calls
            do j = a % row(c1), a % row(c1 + 1) - 1
              a % val(j) = 0.0
            end do
            b(c1) = eps % n(c1)
            a % val(a % dia(c1)) = 1.0
          else
            eps % n(c2) = eps_int
          end if  ! y_plus < 4
        end if    ! rough_walls
      end if      ! wall or wall_flux
    end if        ! c2 < 0
  end do  

  end subroutine
