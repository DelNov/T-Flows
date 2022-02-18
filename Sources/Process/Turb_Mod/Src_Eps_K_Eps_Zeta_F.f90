!==============================================================================!
  subroutine Turb_Mod_Src_Eps_K_Eps_Zeta_F(turb, Nat)
!------------------------------------------------------------------------------!
!   Calculates source terms in equation of dissipation of turbulent energy     !
!   and imposes boundary condition                                             !
!------------------------------------------------------------------------------!
  implicit none
!--------------------------------[Arguments]-----------------------------------!
  type(Turb_Type),   target :: turb
  type(Native_Type), target :: Nat
!---------------------------------[Calling]------------------------------------!
  real :: Roughness_Coefficient
  real :: Tau_Wall_Low_Re
  real :: Y_Plus_Low_Re
  real :: Y_Plus_Rough_Walls
  real :: Tau_Wall_Rough_Walls
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type),  pointer :: Flow
  type(Grid_Type),   pointer :: Grid
  type(Var_Type),    pointer :: u, v, w
  type(Var_Type),    pointer :: kin, eps, zeta, f22, ut, vt, wt
  type(Matrix_Type), pointer :: A
  real,              pointer :: b(:)
  integer                    :: c, s, c1, c2, j
  real                       :: u_tan, u_tau
  real                       :: e_sor, c_11e, ebf
  real                       :: eps_wf, eps_int
  real                       :: fa, u_tau_new, kin_vis, p_kin_int, p_kin_wf
  real                       :: z_o
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
!   left hand s.  A        [kg/s]      | right hand s.   b         [kg*m^2/s^4]!
!------------------------------------------------------------------------------!
!   p_kin = 2*vis_t / density S_ij S_ij                                        !
!   shear = sqrt(2 S_ij S_ij)                                                  !
!------------------------------------------------------------------------------!

  ! Take aliases
  Flow => turb % pnt_flow
  Grid => Flow % pnt_grid
  call Flow % Alias_Momentum(u, v, w)
  call Turb_Mod_Alias_K_Eps_Zeta_F(turb, kin, eps, zeta, f22)
  call Turb_Mod_Alias_Heat_Fluxes (turb, ut, vt, wt)
  call Nat % Alias_Native         (A, b)

  call Time_And_Length_Scale(Grid, turb)

  do c = 1, Grid % n_cells 
    e_sor = Grid % vol(c)/(turb % t_scale(c)+TINY)
    c_11e = c_1e*(1.0 + alpha * ( 1.0/(zeta % n(c)+TINY) ))
    b(c) = b(c) + c_11e * e_sor * turb % p_kin(c)

    ! Fill in a diagonal of coefficient matrix
    A % val(A % dia(c)) =  A % val(A % dia(c))  &
                        + c_2e * e_sor * Flow % density(c)

    ! Add buoyancy (linearly split) to eps equation as required in the t2 model
    if(Flow % buoyancy .eq. THERMALLY_DRIVEN) then
      b(c) = b(c) + max(0.0, c_11e * e_sor * turb % g_buoy(c))
      A % val(A % dia(c)) = A % val(A % dia(c))                         &
                          + max(0.0, -c_11e * e_sor * turb % g_buoy(c)  &
                          / (eps % n(c) + TINY))
    end if
  end do

  !-------------------------------------------------------!
  !   Following block shows density dependent behaviour   !
  !-------------------------------------------------------!

  ! Imposing a boundary condition on wall for eps
  do s = 1, Grid % n_faces
    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)
    if(c2 < 0) then
      kin_vis = Flow % viscosity(c1) / Flow % density(c1)
      if( Grid % Bnd_Cond_Type(c2) .eq. WALL .or.  &
          Grid % Bnd_Cond_Type(c2) .eq. WALLFL) then

        ! Compute tangential velocity component
        u_tan = Flow % U_Tan(s)

        u_tau = c_mu25 * sqrt(kin % n(c1))

        turb % y_plus(c1) = Y_Plus_Low_Re(turb,                  &
                                          u_tau,                 &
                                          Grid % wall_dist(c1),  &
                                          kin_vis)

        turb % tau_wall(c1) = Tau_Wall_Low_Re(turb,               &
                                              Flow % density(c1), &
                                              u_tau,              &
                                              u_tan,              &
                                              turb % y_plus(c1))

        u_tau_new = sqrt(turb % tau_wall(c1) / Flow % density(c1))

        turb % y_plus(c1) = Y_Plus_Low_Re(turb,                  &
                                          u_tau_new,             &
                                          Grid % wall_dist(c1),  &
                                          kin_vis)

        eps_int = 2.0* kin_vis * kin % n(c1)  &
                / Grid % wall_dist(c1)**2
        eps_wf  = c_mu75 * kin % n(c1)**1.5            &
                / (Grid % wall_dist(c1) * kappa)

        ebf = Turb_Mod_Ebf_Momentum(turb, c1)

        p_kin_wf  = turb % tau_wall(c1) * c_mu25 * sqrt(kin % n(c1))  &
                / (Grid % wall_dist(c1) * kappa)

        p_kin_int = turb % vis_t(c1) * Flow % shear(c1)**2

        turb % p_kin(c1) = p_kin_int * exp(-1.0 * ebf) + p_kin_wf  &
                           * exp(-1.0 / ebf)

        fa = min( p_kin_wf * exp(-1.0 / ebf) / (turb % p_kin(c1) + TINY), 1.0)

        eps % n(c1) = (1.0 - fa)**0.5 * eps_int + fa**0.5 * eps_wf

        if(turb % rough_walls) then
          z_o = Roughness_Coefficient(turb, turb % z_o_f(c1))
          z_o = max(Grid % wall_dist(c1)  &
              / (e_log * max(turb % y_plus(c1), 1.0)), z_o)
 
          turb % y_plus(c1) = Y_Plus_Rough_Walls(turb,                  &
                                                 u_tau,                 &
                                                 Grid % wall_dist(c1),  &
                                                 kin_vis)

          turb % tau_wall(c1) = Tau_Wall_Rough_Walls(turb,                  &
                                                     Flow % density(c1),    &
                                                     u_tau,                 &
                                                     u_tan,                 &
                                                     Grid % wall_dist(c1),  &
                                                     z_o)

          p_kin_wf = turb % tau_wall(c1) * c_mu25 * sqrt(kin % n(c1))  &
                   / (kappa * (Grid % wall_dist(c1) + z_o))

          eps_wf = c_mu75 * kin % n(c1)**1.5  &
                 / ((Grid % wall_dist(c1) + z_o) * kappa)

          ebf = Turb_Mod_Ebf_Momentum(turb, c1)

          p_kin_int = turb % vis_t(c1) * Flow % shear(c1)**2

          turb % p_kin(c1) = p_kin_int * exp(-1.0 * ebf)  &
                           + p_kin_wf * exp(-1.0 / ebf)

          fa = min( p_kin_wf * exp(-1.0 / ebf) / (turb % p_kin(c1) + TINY), 1.0)

          eps % n(c1) = (1.0 - fa)**0.5 * eps_int  &
                      + fa**0.5 * eps_wf
        end if  ! rough_walls

        if(turb % y_plus(c1) > 3) then
          ! Adjusting coefficient to fix eps value in near wall calls
          do j = A % row(c1), A % row(c1 + 1) - 1
            A % val(j) = 0.0
          end do
          b(c1) = eps % n(c1)
          A % val(A % dia(c1)) = 1.0
        else
          eps % n(c2) = 2.0* kin_vis * kin % n(c1)  &
                      / Grid % wall_dist(c1)**2
        end if  ! y_plus(c1) < 3
      end if      ! wall or wall_flux
    end if        ! c2 < 0
  end do

  end subroutine
