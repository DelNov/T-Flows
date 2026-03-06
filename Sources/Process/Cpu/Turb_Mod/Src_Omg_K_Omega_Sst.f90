!==============================================================================!
  subroutine Src_Omg_K_Omega_Sst(Turb, Sol)
!------------------------------------------------------------------------------!
!   Computes the source terms in the omega transport equation (k-omega SST),
!   keeping the existing wall-function machinery (tau_wall, y_plus, kappa, etc.)
!   for compatibility with wall functions and heat transfer.
!
!   Omega equation (simplified, SST-style):
!
!     d(rho*omega)/dt + div(rho*u*omega) = div( (mu + mu_t/sigma_w) grad omega)
!                                        + gamma * rho * Pk / mu_t
!                                        - beta  * rho * omega^2
!
!   Implicit linearization of -beta*rho*omega^2:
!     A_dia += beta * rho * omega_old * V
!
!   Wall-function omega (when wall-functions are used):
!     omega_wf = u_tau / (kappa * (y + z0) * sqrt(beta_star))
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Turb_Type),  target :: Turb
  type(Solver_Type), target :: Sol
!------------------------------[Local parameters]------------------------------!
! real, parameter :: BETA1  = 0.0750
! real, parameter :: BETA2  = 0.0828
! real, parameter :: GAMMA1 = 5.0/9.0
! real, parameter :: GAMMA2 = 0.44
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type),  pointer :: Flow
  type(Grid_Type),   pointer :: Grid
  type(Var_Type),    pointer :: u, v, w
  type(Var_Type),    pointer :: kin, omega
  type(Matrix_Type), pointer :: A
  real,              pointer :: b(:)
  integer                    :: s, c, c1, c2, j, reg
  real                       :: u_tan, u_tau
  real                       :: kin_vis
  real                       :: omg_wf, omg_int
  real                       :: p_kin_int, p_kin_wf, ebf, z_o, fa
  real                       :: mu_t
  real                       :: f1, beta_c, gamma_c, dia_coef_tmp
!==============================================================================!

  ! Take aliases
  Flow => Turb % pnt_flow
  Grid => Flow % pnt_grid

  call Flow % Alias_Momentum(u, v, w)

  kin   => Turb % kin
  omega => Turb % omega

  call Sol % Alias_Native(A, b)

  !---------------------------------------------!
  !                                             !
  !   Domain source terms for omega equation    !
  !                                             !
  !---------------------------------------------!
  do c = Cells_In_Domain()

    kin_vis = Flow % viscosity(c) / Flow % density(c)
    mu_t    = max(Turb % vis_t(c), TINY)

! ---- SAFE f1 fetch (po ćeliji) ----
!    if (allocated(Turb%sst_f1)) then
!      if (size(Turb%sst_f1) >= c) then
!        f1 = Turb%sst_f1(c)
!      else
!        f1 = 1.0
!      end if
!    else
!      f1 = 1.0
!    end if

    ! SST blending (f1): use cell-wise beta/gamma
    f1      = Turb % sst_f1(c)
    beta_c  = f1 * Turb % beta1  + (1.0-f1) * Turb % beta2
    gamma_c = f1 * Turb % gamma1 + (1.0-f1) * Turb % gamma2

    ! Positive contribution: gamma * rho * Pk / mu_t
    b(c) = b(c) + gamma_c * Flow % density(c) * Turb % p_kin(c)  &
                 / mu_t * Grid % vol(c)

    ! Negative contribution: -beta * rho * omega^2  (implicit)
    A % val(A % dia(c)) = A % val(A % dia(c))  &
                        + Flow % density(c) * beta_c * omega % n(c)  &
                        * Grid % vol(c)
  end do

  !--------------------------------------------------!
  !                                                  !
  !   Wall function treatment (keep existing logic)  !
  !                                                  !
  !--------------------------------------------------!
  do reg = Boundary_Regions()

    if(Grid % region % type(reg) .eq. WALL .or.  &
       Grid % region % type(reg) .eq. WALLFL) then

      do s = Faces_In_Region(reg)
        c1 = Grid % faces_c(1,s)
        c2 = Grid % faces_c(2,s)

        Assert(c2 < 0)  ! just to make sure

        ! Kinematic viscosity at the near-wall cell
        kin_vis = Flow % viscosity(c1) / Flow % density(c1)

        ! Set up roughness coefficient (keep your existing function)
        z_o = Turb % Roughness_Coeff(c1, c2)

        ! Tangential velocity at the wall face
        u_tan = Flow % U_Tan(s)

        ! Friction velocity from k (same as in eps wall-function routine)
        u_tau = Turb % c_mu25 * sqrt(max(kin % n(c1), 0.0))

        ! y+ (keep your existing rough-wall function)
        Turb % y_plus(c1) = Turb % Y_Plus_Rough_Walls(    &
                                   u_tau,                 &
                                   Grid % wall_dist(c1),  &
                                   kin_vis,               &
                                   z_o)

        ! Wall-function omega and interior omega; blend like eps routine
        omg_wf = u_tau / ((Grid % wall_dist(c1) + z_o)    &
                 * Turb % kappa * sqrt(Turb % beta_star))


        omg_int = 60.0 * kin_vis / (Turb % beta1                   &
                  * max(Grid % wall_dist(c1),TINY)**2)

        ! If you use wall-functions, enforce boundary value on the wall face

        ebf = Turb % Ebf_Momentum(c1)

        p_kin_wf = Turb % tau_wall(c1) * Turb % c_mu25 * sqrt(kin % n(c1))  &
                 / ((Grid % wall_dist(c1) + z_o) * Turb % kappa)

        p_kin_int = Turb % vis_t(c1) * Flow % shear(c1)**2

        Turb % p_kin(c1) = p_kin_int * exp(-1.0 * ebf) + p_kin_wf  &
                         * exp(-1.0 / ebf)

        fa = min(p_kin_wf * exp(-1.0 / ebf) / (Turb % p_kin(c1) + TINY), 1.0)


        if(Turb % y_plus(c1) > 3) then
          omega % n(c1) = omg_wf !(1.0 - fa)**0.5 * omg_int + fa**0.5 * omg_wf
          omega % n(c2) = 0.0

          dia_coef_tmp = A % val(A % dia(c1))

          ! Adjusting coefficient to fix eps value in near wall calls
          do j = A % row(c1), A % row(c1 + 1) - 1
            A % val(j) = 0.0
          end do

          b(c1) = omega % n(c1) * dia_coef_tmp
          A % val(A % dia(c1)) = dia_coef_tmp
        else

          omega % n(c2) = 60.0 * kin_vis     &
                        / (Turb % beta1 * max(Grid % wall_dist(c1),TINY)**2)
        end if  ! y_plus(c1) < 3

      end do

    end if

  end do

  ! Ensure omega stays positive in domain + buffers
  do c = Cells_In_Domain_And_Buffers()
    omega % n(c) = max(omega % n(c), TINY)
  end do

  end subroutine Src_Omg_K_Omega_Sst
!==============================================================================!
