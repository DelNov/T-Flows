!==============================================================================!
  subroutine Src_Eps_K_Eps(Turb, Sol)
!------------------------------------------------------------------------------!
!   Computes the source terms in the eps transport equation,                   !
!   wall shear stress (wall function approuch)                                 !
!------------------------------------------------------------------------------!
!
!  The form of the term being discretized:
!                                                                              !
!     /                                                                        !
!    |                                                                         !
!    | ( density (c_1e eps/kin Gk - c_2e eps^2/kin) ) dV                       !
!    |                                                                         !
!   /                                                                          !
!                                                                              !
!   assigns epsilon from the wall function:                                    !
!                                                                              !
!   Eps_w = Cmu^(3/4)* Kin^(3/2)/(Kappa/y)                                     !
!                                                                              !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Turb_Type),  target :: Turb
  type(Solver_Type), target :: Sol
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type),  pointer :: Flow
  type(Grid_Type),   pointer :: Grid
  type(Var_Type),    pointer :: u, v, w
  type(Var_Type),    pointer :: kin, eps
  type(Matrix_Type), pointer :: A
  real,              pointer :: b(:)
  integer                    :: s, c, c1, c2, j, reg
  real                       :: u_tan, u_tau
  real                       :: re_t, f_mu, fa, kin_vis
  real                       :: eps_wf, eps_int, y_star, dia_coef_tmp
  real                       :: p_kin_int, p_kin_wf, ebf, z_o
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
  call Turb % Alias_K_Eps   (kin, eps)
  call Sol % Alias_Native   (A, b)

  do c = Cells_In_Domain()
    kin_vis =  Flow % viscosity(c) / Flow % density(c)

    ! Positive contribution:
    b(c) = b(c)  &
         + Turb % c_1e * Turb % p_kin(c) * eps % n(c)  &
         / kin % n(c) * Grid % vol(c)

    ! Negative contribution:
    re_t = kin % n(c)*kin % n(c)/(kin_vis*eps % n(c))
    y_star = sqrt(sqrt(kin_vis * eps % n(c))) *     &
             Grid % wall_dist(c)/kin_vis
    f_mu = (1.0 - exp(-y_star/3.1))**2              &
         * (1.0 - 0.3*exp(-(re_t/6.5)*(re_t/6.5)))

    f_mu = min(f_mu,1.0)

    A % val(A % dia(c)) = A % val(A % dia(c))                             &
                +    Flow % density(c) * f_mu * Turb % c_2e * eps % n(c)  &
                   / kin % n(c) * Grid % vol(c)

    ! Buoyancy contribution
    if(Flow % buoyancy .eq. THERMALLY_DRIVEN) then
      b(c) = b(c) + max(0.0, Turb % c_1e * Turb % g_buoy(c) &
                    * eps % n(c) / kin % n(c) * Grid % vol(c))
      A % val(A % dia(c)) = A % val(A % dia(c))                &
                          + max(0.0,(-Turb % c_1e * Turb % g_buoy(c)  &
                          * eps % n(c)                         &
                          / kin % n(c) * Grid % vol(c))        &
                          / (eps % n(c) + TINY))
    end if

  end do

  ! Imposing a boundary condition on wall for eps
  do reg = Boundary_Regions()
    if(Grid % region % type(reg) .eq. WALL .or.  &
       Grid % region % type(reg) .eq. WALLFL) then
      do s = Faces_In_Region(reg)
        c1 = Grid % faces_c(1,s)
        c2 = Grid % faces_c(2,s)

        Assert(c2 < 0)  ! just to make sure

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
