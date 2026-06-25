!==============================================================================!
  subroutine Src_Eps_K_Eps_Zeta_F(Turb, Sol)
!------------------------------------------------------------------------------!
!   Calculates source terms in equation of dissipation of turbulent energy     !
!   and imposes boundary condition                                             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Turb_Type),  target :: Turb
  type(Solver_Type), target :: Sol
!------------------------------[Local parameters]------------------------------!
  real, parameter :: RE_T_CUT  =    10.0
  real, parameter :: RE_T_1    =   600.0
  real, parameter :: RE_T_4    = 20000.0
  real, parameter :: ALPHA_MIN =     0.012
  real, parameter :: ALPHA_MAX =     0.018
  real, parameter :: CUB_A3    =    -1.24081684e-16
  real, parameter :: CUB_A2    =     4.73275447e-12
  real, parameter :: CUB_A1    =     2.62949932e-7
  real, parameter :: CUB_A0    =     1.18405531e-2
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type),  pointer :: Flow
  type(Grid_Type),   pointer :: Grid
  type(Var_Type),    pointer :: u, v, w
  type(Var_Type),    pointer :: kin, eps, zeta, f22, ut, vt, wt
  type(Matrix_Type), pointer :: A
  real,              pointer :: b(:)
  integer                    :: c, s, c1, c2, j, reg
  real                       :: u_tan, u_tau
  real                       :: e_sor, c_11e, ebf
  real                       :: eps_wf, eps_int
  real                       :: fa, kin_vis, p_kin_int, p_kin_wf
  real                       :: z_o, dia_coef_tmp
  real                       :: re_t_loc, re_t_avg, sum_re_t, sum_vol
  real                       :: alpha_zeta, l_sgs_d, l_rans_d, alpha_d
!------------------------------------------------------------------------------!
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
!==============================================================================!

  ! Take aliases
  Flow => Turb % pnt_flow
  Grid => Flow % pnt_grid
  call Flow % Alias_Momentum(u, v, w)
  call Turb % Alias_K_Eps_Zeta_F(kin, eps, zeta, f22)
  call Turb % Alias_Heat_Fluxes (ut, vt, wt)
  call Sol % Alias_Native       (A, b)

  call Turb % Time_And_Length_Scale(Grid)

  !------------------------------------------------------------------!
  !   Reynolds-number-dependent coefficient in C_e1 correction       !
  !                                                                  !
  !   C_e1 = c_1e * (1 + alpha_zeta / zeta)                          !
  !                                                                  !
  !   alpha_zeta is computed once for the whole domain. First,       !
  !   the maximum wall distance is found. Then Re_t is averaged      !
  !   in the layer y <= 0.1*y_max, excluding cells with Re_t <= 10.  !
  !   Finally, alpha_zeta is obtained from a cubic interpolation     !
  !   through the points:                                            !
  !      (600,   0.012)                                              !
  !      (4135,  0.013)                                              !
  !      (9000,  0.014)                                              !
  !      (20000, 0.018)                                              !
  !------------------------------------------------------------------!
  sum_re_t = 0.0
  sum_vol  = 0.0

  do c = Cells_In_Domain()

    kin_vis = Flow % viscosity(c) / Flow % density(c)

    ! Local turbulent Reynolds number:
    !
    !   Re_t = k^2 / (nu * eps)
    !
    re_t_loc = max(kin % n(c), 0.0)**2                         &
             / (kin_vis * max(eps % n(c), TINY) + TINY)

    if(Turb % model == HYBRID_LES_RANS) then
      l_sgs_d  = Grid % vol(c)**ONE_THIRD
      l_rans_d = Turb % kappa * Grid % wall_dist(c)
      alpha_d  = l_rans_d/l_sgs_d

      if( alpha_d < Turb % c_hyb ) then
        if(re_t_loc > RE_T_CUT) then
          sum_re_t = sum_re_t + Grid % vol(c) * re_t_loc
          sum_vol  = sum_vol  + Grid % vol(c)
        end if
      end if
    else
      if(re_t_loc > RE_T_CUT) then
        sum_re_t = sum_re_t + Grid % vol(c) * re_t_loc
        sum_vol  = sum_vol  + Grid % vol(c)
      end if
    end if

  end do

  call Global % Sum_Real(sum_re_t)
  call Global % Sum_Real(sum_vol)

  if(sum_vol > TINY) then
    re_t_avg = sum_re_t / sum_vol
  else
    re_t_avg = RE_T_1
  end if

  if(re_t_avg <= RE_T_1) then
    alpha_zeta = ALPHA_MIN
  else if(re_t_avg >= RE_T_4) then
    alpha_zeta = ALPHA_MAX
  else
    alpha_zeta = (  (CUB_A3 * re_t_avg + CUB_A2) * re_t_avg  &
                  + CUB_A1) * re_t_avg + CUB_A0
    alpha_zeta = max(ALPHA_MIN, min(ALPHA_MAX, alpha_zeta))
  end if

  do c = Cells_In_Domain()
    kin_vis = Flow % viscosity(c) / Flow % density(c)
    e_sor = Grid % vol(c) / (Turb % t_scale(c)+TINY)
    c_11e = Turb % c_1e*(1.0 + alpha_zeta * (1.0/(zeta % n(c)+TINY)))
    b(c) = b(c) + c_11e * e_sor * Turb % p_kin(c)

    ! Fill in a diagonal of coefficient matrix
    A % val(A % dia(c)) = A % val(A % dia(c))  &
                        + Turb % c_2e * e_sor * Flow % density(c)

    ! Add buoyancy (linearly split) to eps equation as required in the t2 model
    if(Flow % buoyancy .eq. THERMALLY_DRIVEN) then
      b(c) = b(c) + max(0.0, c_11e * e_sor * Turb % g_buoy(c))
      A % val(A % dia(c)) = A % val(A % dia(c))                         &
                          + max(0.0, -c_11e * e_sor * Turb % g_buoy(c)  &
                          / (eps % n(c) + TINY))
    end if
  end do

  !-------------------------------------------------------!
  !   Following block shows density dependent behaviour   !
  !-------------------------------------------------------!

  ! Imposing a boundary condition on wall for eps
  do reg = Boundary_Regions()
    if(Grid % region % type(reg) .eq. WALL .or.  &
       Grid % region % type(reg) .eq. WALLFL) then
      do s = Faces_In_Region(reg)
        c1 = Grid % faces_c(1,s)
        c2 = Grid % faces_c(2,s)

        Assert(c2 < 0)

        kin_vis = Flow % viscosity(c1) / Flow % density(c1)

        ! Set up roughness coefficient
        z_o = Turb % Roughness_Coeff(c1, c2)

        ! Compute tangential velocity component
        u_tan = Flow % U_Tan(s)

        u_tau = Turb % c_mu25 * sqrt(kin % n(c1))

        Turb % y_plus(c1) = Turb % Y_Plus_Rough_Walls(    &
                                   u_tau,                 &
                                   Grid % wall_dist(c1),  &
                                   kin_vis,               &
                                   z_o)

        eps_int = 2.0 * kin_vis * kin % n(c1)  &
                / Grid % wall_dist(c1)**2

        eps_wf  = Turb % c_mu75 * kin % n(c1)**1.5   &
                / ((Grid % wall_dist(c1) + z_o) * Turb % kappa)

        ebf = Turb % Ebf_Momentum(c1)

        p_kin_wf = Turb % tau_wall(c1) * Turb % c_mu25 * sqrt(kin % n(c1))  &
                 / ((Grid % wall_dist(c1) + z_o) * Turb % kappa)

        p_kin_int = Turb % vis_t(c1) * Flow % shear(c1)**2

        Turb % p_kin(c1) = p_kin_int * exp(-1.0 * ebf) + p_kin_wf  &
                         * exp(-1.0 / ebf)

        fa = min(p_kin_wf * exp(-1.0 / ebf) / (Turb % p_kin(c1) + TINY), 1.0)

        eps % n(c1) = (1.0 - fa)**0.5 * eps_int + fa**0.5 * eps_wf

        if(Turb % y_plus(c1) > 3) then

          dia_coef_tmp = A % val(A % dia(c1))

          ! Adjusting coefficient to fix eps value in near wall calls
          do j = A % row(c1), A % row(c1 + 1) - 1
            A % val(j) = 0.0
          end do

          b(c1) = eps % n(c1) * dia_coef_tmp
          A % val(A % dia(c1)) = dia_coef_tmp

        else

          eps % n(c2) = 2.0 * kin_vis * kin % n(c1)  &
                      / Grid % wall_dist(c1)**2
        end if  ! y_plus(c1) < 3
      end do    ! faces in regions
    end if      ! region is WALL or WALLFL
  end do        ! through regions

  end subroutine
