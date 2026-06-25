!==============================================================================!
  subroutine User_Mod_Save_Results(Flow, Turb, Vof, Swarm, domain)
!------------------------------------------------------------------------------!
!   Averages results in homogeneous (z) planes and writes them to a .dat file. !
!   Works for all turbulence models: RANS, LES, and hybrid.                   !
!                                                                              !
!   Statistics handling:                                                       !
!     Turb % statistics == .true.  -> LES/HYB: use time-averaged means and    !
!                                     resolved Reynolds stresses               !
!     Turb % statistics == .false. -> RANS (or LES/HYB before stat window):   !
!                                     use instantaneous field                  !
!                                                                              !
!   Mean arrays (u_mean, uu_res, …) are ONLY allocated when                   !
!   Turb % statistics == .true. (see Create_Turb.f90), therefore use_stats is  !
!   additionally guarded by allocated(...) checks before accessing them.        !
!                                                                              !
!   Friction-coefficient reference: Zanoun et al. implicit correlation         !
!     1/sqrt(Cf) = 1.911 * ln(Re * sqrt(Cf)) – 1.282                         !
!   Solved by Newton iteration (Dean used as initial guess only).              !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type), target :: Flow
  type(Turb_Type),  target :: Turb
  type(Vof_Type),   target :: Vof
  type(Swarm_Type), target :: Swarm
  integer, optional        :: domain
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
  type(Bulk_Type), pointer :: bulk
  type(Var_Type),  pointer :: u, v, w, kin, eps, zeta, f22, omega, t
  integer                  :: n_prob, i, c, s, c1, c2, n_points, fu, iter, pl
  character(SL)            :: res_name
  !---------------------------------------------------------------------------!
  !   Local averaging arrays (always allocated, selectively filled)           !
  !---------------------------------------------------------------------------!
  real, allocatable :: wall_p(:)                 ! mean wall distance
  real, allocatable :: u_p(:), v_p(:), w_p(:)   ! mean velocity components
  real, allocatable :: y_plus_p(:)               ! mean y+
  real, allocatable :: uu_p(:), vv_p(:), ww_p(:) ! normal resolved stresses
  real, allocatable :: uw_p(:)                   ! shear: resolved (stats) or
                                                  ! modelled (RANS)
  real, allocatable :: uw_mod_p(:)               ! modelled shear (LES/HYB)
  real, allocatable :: kin_p(:), eps_p(:)        ! k, eps
  real, allocatable :: omg_p(:)                  ! omega (K_OMEGA_SST)
  real, allocatable :: f22_p(:), zeta_p(:)       ! f22, Re_t (K_EPS_ZETA_F)
  real, allocatable :: vis_t_p(:)                ! eddy viscosity / visc_const
  real, allocatable :: t_p(:)                    ! temperature
  real, allocatable :: t2_p(:)                   ! temperature variance
  real, allocatable :: ut_p(:), vt_p(:), wt_p(:) ! turbulent heat fluxes
  integer, allocatable :: n_count(:)             ! cell count per plane
  !---------------------------------------------------------------------------!
  real    :: ubulk, re, cf, cf_zanoun, error, pr
  real    :: zan_fun, zan_der
  real    :: u_tau_p, t_tau, t_wall, nu_mean
  real    :: dens_const, visc_const, capa_const, cond_const
  real    :: kin_vis, u_vel, v_vel, w_vel
  logical :: use_stats, use_heat_stats
!==============================================================================!

  !------------------------------------------------------------------!
  !   Return at t=0: the solution hasn't developed yet               !
  !------------------------------------------------------------------!
  if(Time % Curr_Dt() .eq. 0) return

  !-----------------------------!
  !   Aliases and properties    !
  !-----------------------------!
  Grid  => Flow % pnt_grid
  bulk  => Flow % bulk
  call Flow % Alias_Momentum(u, v, w)
  call Turb % Alias_K_Eps_Zeta_F(kin, eps, zeta, f22)
  call Flow % Alias_Energy(t)
  omega => Turb % omega

  call Control % Mass_Density        (dens_const)
  call Control % Dynamic_Viscosity   (visc_const)
  call Control % Heat_Capacity       (capa_const)
  call Control % Thermal_Conductivity(cond_const)

  !------------------------------------------------------------------!
  !   Decide whether to use time-averaged statistics.                !
  !   The statistics flag alone is not enough: mean / resolved       !
  !   arrays are accessed only if they are really allocated.         !
  !   If not, the routine safely falls back to instantaneous fields. !
  !------------------------------------------------------------------!
  use_stats = Turb % statistics             .and.  &
              allocated(Turb % u_mean)      .and.  &
              allocated(Turb % v_mean)      .and.  &
              allocated(Turb % w_mean)      .and.  &
              allocated(Turb % uu_res)      .and.  &
              allocated(Turb % vv_res)      .and.  &
              allocated(Turb % ww_res)      .and.  &
              allocated(Turb % uw_res)

  use_heat_stats = use_stats
  if(Flow % heat_transfer) then
    use_heat_stats = use_stats              .and.  &
                     allocated(Turb % t_mean)      .and.  &
                     allocated(Turb % t2_res)      .and.  &
                     allocated(Turb % ut_res)      .and.  &
                     allocated(Turb % vt_res)      .and.  &
                     allocated(Turb % wt_res)
  end if

  ubulk    = bulk % flux_x / bulk % area_x
  t_wall   = 0.0
  nu_mean  = 0.0
  n_points = 0

  call File % Set_Name(res_name,                     &
                       time_step = Time % Curr_Dt(), &
                       extension = '-res.dat')

  ! Number of averaging intervals between consecutive z-coordinate planes
  n_prob = Grid % n_z_planes - 1

  !-------------------------------------------------------------------------!
  !   Allocate local arrays (all of them; unused ones simply stay at zero)   !
  !-------------------------------------------------------------------------!
  allocate(n_count  (n_prob));  n_count   = 0
  allocate(wall_p   (n_prob));  wall_p    = 0.0
  allocate(u_p      (n_prob));  u_p       = 0.0
  allocate(v_p      (n_prob));  v_p       = 0.0
  allocate(w_p      (n_prob));  w_p       = 0.0
  allocate(y_plus_p (n_prob));  y_plus_p  = 0.0
  allocate(uu_p     (n_prob));  uu_p      = 0.0
  allocate(vv_p     (n_prob));  vv_p      = 0.0
  allocate(ww_p     (n_prob));  ww_p      = 0.0
  allocate(uw_p     (n_prob));  uw_p      = 0.0
  allocate(uw_mod_p (n_prob));  uw_mod_p  = 0.0
  allocate(kin_p    (n_prob));  kin_p     = 0.0
  allocate(eps_p    (n_prob));  eps_p     = 0.0
  allocate(omg_p    (n_prob));  omg_p     = 0.0
  allocate(f22_p    (n_prob));  f22_p     = 0.0
  allocate(zeta_p   (n_prob));  zeta_p    = 0.0
  allocate(vis_t_p  (n_prob));  vis_t_p   = 0.0

  if(Flow % heat_transfer) then
    allocate(t_p (n_prob));  t_p  = 0.0
    allocate(t2_p(n_prob));  t2_p = 0.0
    allocate(ut_p(n_prob));  ut_p = 0.0
    allocate(vt_p(n_prob));  vt_p = 0.0
    allocate(wt_p(n_prob));  wt_p = 0.0
  end if

  !---------------------------------------------------!
  !   Accumulate quantities in each z-plane interval  !
  !---------------------------------------------------!
  do i = 1, n_prob
    do c = Cells_In_Domain()
      if(Grid % zc(c) > Grid % z_coord_plane(i)   .and.  &
         Grid % zc(c) < Grid % z_coord_plane(i+1)) then

        kin_vis = Flow % viscosity(c) / Flow % density(c)

        !--- Choose mean or instantaneous velocity ---------------------------!
        if(use_stats) then
          u_vel = Turb % u_mean(c)
          v_vel = Turb % v_mean(c)
          w_vel = Turb % w_mean(c)
        else
          u_vel = u % n(c)
          v_vel = v % n(c)
          w_vel = w % n(c)
        end if

        wall_p  (i) = wall_p  (i) + Grid % wall_dist(c)
        u_p     (i) = u_p     (i) + u_vel
        v_p     (i) = v_p     (i) + v_vel
        w_p     (i) = w_p     (i) + w_vel
        y_plus_p(i) = y_plus_p(i) + Turb % y_plus(c)

        !=====================================================================!
        !   LES / HYB branch: use time-averaged statistics                    !
        !=====================================================================!
        if(use_stats) then

          !--- Resolved Reynolds stresses (fluctuation = res - mean²) --------!
          uu_p(i) = uu_p(i) + Turb % uu_res(c) - u_vel * u_vel
          vv_p(i) = vv_p(i) + Turb % vv_res(c) - v_vel * v_vel
          ww_p(i) = ww_p(i) + Turb % ww_res(c) - w_vel * w_vel
          uw_p(i) = uw_p(i) + Turb % uw_res(c) - u_vel * w_vel

          !--- Model-specific modelled stress and eddy viscosity -------------!
          if(Turb % model == HYBRID_LES_RANS) then
            ! RANS part carries k and eps; LES part provides vis_t_sgs
            kin_p   (i) = kin_p   (i) + Turb % kin_mean(c)
            eps_p   (i) = eps_p   (i) + Turb % eps_mean(c)
            vis_t_p (i) = vis_t_p (i) + Turb % vis_t_eff(c) / visc_const
            uw_mod_p(i) = uw_mod_p(i) + Turb % vis_t_eff(c)  &
                                       * (u % z(c) + w % x(c))

          else if(Turb % model == DES_SPALART) then
            vis_t_p (i) = vis_t_p (i) + Turb % vis_t(c) / visc_const
            uw_mod_p(i) = uw_mod_p(i) + Turb % vis_t(c)  &
                                       * (u % z(c) + w % x(c))

          else if(Turb % model /= DNS              .and.  &
                  Turb % model /= NO_TURBULENCE_MODEL) then
            ! All other LES sub-grid models that carry vis_t
            vis_t_p (i) = vis_t_p (i) + Turb % vis_t(c) / visc_const
            uw_mod_p(i) = uw_mod_p(i) + Turb % vis_t(c)  &
                                       * (u % z(c) + w % x(c))
          end if

          !--- Heat transfer ---------------------------------------------------!
          if(Flow % heat_transfer) then
            if(use_heat_stats) then
              t_p (i) = t_p (i) + Turb % t_mean(c)
              ! t2_res = <T²>, so variance = t2_res - <T>²
              ! (t2_res is allocated for all statistics-enabled models with heat
              !  transfer; t2_mean is only available for HYBRID_LES_RANS / RANS)
              t2_p(i) = t2_p(i) + Turb % t2_res(c)         &
                                - Turb % t_mean(c) * Turb % t_mean(c)
              ut_p(i) = ut_p(i) + Turb % ut_res(c)         &
                                - u_vel * Turb % t_mean(c)
              vt_p(i) = vt_p(i) + Turb % vt_res(c)         &
                                - v_vel * Turb % t_mean(c)
              wt_p(i) = wt_p(i) + Turb % wt_res(c)         &
                                - w_vel * Turb % t_mean(c)
            else
              t_p(i) = t_p(i) + t % n(c)
            end if
          end if

        !=====================================================================!
        !   RANS branch (or LES/HYB before statistics window opens):          !
        !   use instantaneous field; no resolved stresses available           !
        !=====================================================================!
        else

          !--- Modelled shear stress (dynamic units) ---------------------------!
          !    Divided by density during non-dimensionalisation below          !
          if(Turb % model == HYBRID_LES_RANS) then
            uw_p(i) = uw_p(i) + Turb % vis_t_eff(c) * (u % z(c) + w % x(c))
          else
            uw_p(i) = uw_p(i) + Turb % vis_t(c) * (u % z(c) + w % x(c))
          end if

          !--- Model-specific quantities --------------------------------------!
          select case(Turb % model)

            case(K_EPS_ZETA_F)
              kin_p  (i) = kin_p  (i) + kin   % n(c)
              eps_p  (i) = eps_p  (i) + eps   % n(c)
              f22_p  (i) = f22_p  (i) + f22   % n(c)
              ! Store turbulent Reynolds number Re_t = k²/(ν·ε) in zeta_p
              if(eps % n(c) > 0.0) then
                zeta_p (i) = zeta_p (i) + kin % n(c) * kin % n(c)  &
                                        / (kin_vis * eps % n(c))
              end if
              vis_t_p(i) = vis_t_p(i) + Turb % vis_t(c) / visc_const

            case(K_EPS)
              kin_p  (i) = kin_p  (i) + kin   % n(c)
              eps_p  (i) = eps_p  (i) + eps   % n(c)
              vis_t_p(i) = vis_t_p(i) + Turb  % vis_t(c) / visc_const

            case(K_OMEGA_SST)
              kin_p  (i) = kin_p  (i) + kin   % n(c)
              omg_p  (i) = omg_p  (i) + omega % n(c)
              vis_t_p(i) = vis_t_p(i) + Turb  % vis_t(c) / visc_const

            case(SPALART_ALLMARAS)
              vis_t_p(i) = vis_t_p(i) + Turb  % vis_t(c) / visc_const

          end select

          if(Flow % heat_transfer) t_p(i) = t_p(i) + t % n(c)

        end if  ! use_stats

        n_count(i) = n_count(i) + 1

      end if
    end do
  end do

  !------------------------------------------!
  !   Global sums across MPI ranks           !
  !------------------------------------------!
  do pl = 1, n_prob
    call Global % Sum_Int (n_count  (pl))
    call Global % Sum_Real(wall_p   (pl))
    call Global % Sum_Real(u_p      (pl))
    call Global % Sum_Real(v_p      (pl))
    call Global % Sum_Real(w_p      (pl))
    call Global % Sum_Real(y_plus_p (pl))
    call Global % Sum_Real(uu_p     (pl))
    call Global % Sum_Real(vv_p     (pl))
    call Global % Sum_Real(ww_p     (pl))
    call Global % Sum_Real(uw_p     (pl))
    call Global % Sum_Real(uw_mod_p (pl))
    call Global % Sum_Real(kin_p    (pl))
    call Global % Sum_Real(eps_p    (pl))
    call Global % Sum_Real(omg_p    (pl))
    call Global % Sum_Real(f22_p    (pl))
    call Global % Sum_Real(zeta_p   (pl))
    call Global % Sum_Real(vis_t_p  (pl))
    if(Flow % heat_transfer) then
      call Global % Sum_Real(t_p  (pl))
      call Global % Sum_Real(t2_p (pl))
      call Global % Sum_Real(ut_p (pl))
      call Global % Sum_Real(vt_p (pl))
      call Global % Sum_Real(wt_p (pl))
    end if
  end do
  call Global % Wait

  !------------------------------------------!
  !   Divide to get per-plane averages       !
  !------------------------------------------!
  do i = 1, n_prob
    if(n_count(i) .ne. 0) then
      wall_p  (i) = wall_p  (i) / n_count(i)
      u_p     (i) = u_p     (i) / n_count(i)
      v_p     (i) = v_p     (i) / n_count(i)
      w_p     (i) = w_p     (i) / n_count(i)
      y_plus_p(i) = y_plus_p(i) / n_count(i)
      uu_p    (i) = uu_p    (i) / n_count(i)
      vv_p    (i) = vv_p    (i) / n_count(i)
      ww_p    (i) = ww_p    (i) / n_count(i)
      uw_p    (i) = uw_p    (i) / n_count(i)
      uw_mod_p(i) = uw_mod_p(i) / n_count(i)
      kin_p   (i) = kin_p   (i) / n_count(i)
      eps_p   (i) = eps_p   (i) / n_count(i)
      omg_p   (i) = omg_p   (i) / n_count(i)
      f22_p   (i) = f22_p   (i) / n_count(i)
      zeta_p  (i) = zeta_p  (i) / n_count(i)
      vis_t_p (i) = vis_t_p (i) / n_count(i)
      if(Flow % heat_transfer) then
        t_p (i) = t_p (i) / n_count(i)
        t2_p(i) = t2_p(i) / n_count(i)
        ut_p(i) = ut_p(i) / n_count(i)
        vt_p(i) = vt_p(i) / n_count(i)
        wt_p(i) = wt_p(i) / n_count(i)
      end if
    end if
  end do

  !------------------------------------------------------------------!
  !   Friction velocity                                              !
  !   y+(1) > 5  → pressure-drop formula (reliable in log region)   !
  !   y+(1) ≤ 5  → viscous-sublayer formula using full |U| near wall !
  !------------------------------------------------------------------!
  if(y_plus_p(1) > 5.0) then
    u_tau_p = sqrt(max(abs(bulk % p_drop_x),  &
                       abs(bulk % p_drop_y),  &
                       abs(bulk % p_drop_z)) / dens_const)
  else
    u_tau_p = sqrt((visc_const                                     &
                   * sqrt(u_p(1)**2 + v_p(1)**2 + w_p(1)**2)      &
                   / wall_p(1)) / dens_const)
  end if

  if(u_tau_p .eq. 0.0) then
    if(First_Proc()) write(*,'(a)')  &
      '# Warning: friction velocity is zero in User_Mod_Save_Results!'
    return
  end if

  !------------------------------------------------------------------!
  !   Zanoun et al. skin-friction correlation for plane channel flow  !
  !     1/sqrt(Cf) = 1.911*ln(Re*sqrt(Cf)) - 1.282                   !
  !   Implicit in Cf → solve by Newton iteration.                     !
  !   Dean's formula ( Cf = 0.073 Re^{-0.25} ) is used as seed only. !
  !------------------------------------------------------------------!
  re        = abs(dens_const * ubulk * 2.0 / visc_const)
  cf        = u_tau_p**2 / (0.5 * ubulk**2)
  cf_zanoun = 0.073 * re**(-0.25)         ! Dean as initial guess

  do iter = 1, 50
    zan_fun = 1.0/sqrt(cf_zanoun)                          &
            - 1.911*log(re*sqrt(cf_zanoun)) + 1.282
    zan_der = -0.5 * cf_zanoun**(-1.5)                     &
              - 0.5 * 1.911 / cf_zanoun
    cf_zanoun = cf_zanoun - zan_fun / zan_der
    if(abs(zan_fun) < 1.0e-8) exit
  end do

  error = abs(cf_zanoun - cf) / cf_zanoun * 100.0

  !------------------------------------------------------------------!
  !   Friction temperature (only with heat transfer)                 !
  !------------------------------------------------------------------!
  if(Flow % heat_transfer) then
    t_tau  = Flow % heat_flux / (dens_const * capa_const * u_tau_p)
    t_wall = 0.0
    do s = 1, Grid % n_faces
      c1 = Grid % faces_c(1,s)
      c2 = Grid % faces_c(2,s)
      if(c2 < 0) then
        if(Grid % Bnd_Cond_Type(c2) .eq. WALL .or.  &
           Grid % Bnd_Cond_Type(c2) .eq. WALLFL) then
          ! Use mean temperature at wall if statistics are collected,
          ! otherwise use the instantaneous boundary value
          if(use_heat_stats) then
            t_wall = t_wall + Turb % t_mean(c2)
          else
            t_wall = t_wall + t % n(c2)
          end if
          n_points = n_points + 1
        end if
      end if
    end do
    call Global % Sum_Real(t_wall)
    call Global % Sum_Int (n_points)
    call Global % Wait
    if(n_points > 0) then
      t_wall = t_wall / n_points
    else
      t_wall = 0.0
    end if
  end if

  !------------------------------------------------------------------!
  !   Open output file and write header                              !
  !------------------------------------------------------------------!
  if(First_Proc()) then

  call File % Open_For_Writing_Ascii(res_name, fu)

  pr = visc_const * capa_const / cond_const

  write(fu,'(a1,(a12,e12.6))') '#', 'Density  = ', dens_const
  write(fu,'(a1,(a12,e12.6))') '#', 'Ubulk    = ', ubulk
  write(fu,'(a1,(a12,e12.6))') '#', 'Re       = ', re
  write(fu,'(a1,(a12,e12.6))') '#', 'Re_tau   = ',  &
                                    dens_const * u_tau_p / visc_const
  write(fu,'(a1,(a12,e12.6))') '#', 'Cf       = ', 2.0*(u_tau_p/ubulk)**2
  write(fu,'(a1,(a12,f12.6))') '#', 'Utau     = ', u_tau_p
  write(fu,'(a1,(a12,f12.6,a2,a26))') '#', 'Cf_error = ', error, ' %',  &
                                      'Zanoun formula is used.'

  !------------------------------------------------------------------!
  !   Non-dimensionalise                                             !
  !                                                                  !
  !   Shear stress convention:                                        !
  !     use_stats  (LES/HYB)  → uw is kinematic (m²/s²),            !
  !                              divide by u_tau²                    !
  !     .not. use_stats (RANS) → uw = vis_t*(du/dz+dw/dx) is        !
  !                              dynamic (Pa), divide by ρ u_tau²   !
  !------------------------------------------------------------------!
  do i = 1, n_prob
    wall_p  (i) = dens_const * wall_p(i) * u_tau_p / visc_const  ! → y+
    u_p     (i) = u_p     (i) / u_tau_p
    v_p     (i) = v_p     (i) / u_tau_p
    w_p     (i) = w_p     (i) / u_tau_p
    kin_p   (i) = kin_p   (i) / u_tau_p**2
    eps_p   (i) = eps_p   (i) * visc_const / (u_tau_p**4 * dens_const)
    f22_p   (i) = f22_p   (i) * visc_const / u_tau_p**2

    if(use_stats) then
      ! Resolved stresses are already kinematic (u'·w' in m²/s²)
      uu_p    (i) = uu_p    (i) / u_tau_p**2
      vv_p    (i) = vv_p    (i) / u_tau_p**2
      ww_p    (i) = ww_p    (i) / u_tau_p**2
      uw_p    (i) = uw_p    (i) / u_tau_p**2
      ! Modelled stress is dynamic (Pa) → divide by ρ as well
      uw_mod_p(i) = uw_mod_p(i) / (u_tau_p**2 * dens_const)
    else
      ! RANS modelled stress is dynamic (Pa)
      uw_p    (i) = uw_p    (i) / (u_tau_p**2 * dens_const)
    end if

    if(Flow % heat_transfer) then
      t_p (i) = (t_wall - t_p(i)) / t_tau
      t2_p(i) = t2_p(i) / (t_tau * t_tau)
      ut_p(i) = ut_p(i) / (u_tau_p * t_tau)
      vt_p(i) = vt_p(i) / (u_tau_p * t_tau)
      wt_p(i) = wt_p(i) / (u_tau_p * t_tau)
    end if
  end do

  !==================================================================!
  !   Write data columns                                             !
  !==================================================================!

  if(use_stats) then

    !----------------------------------------------------------------!
    !   LES / HYB with time-averaged statistics                      !
    !----------------------------------------------------------------!

    if(Turb % model == HYBRID_LES_RANS) then

      !--- Hybrid LES-RANS: resolved + modelled k and shear --------!
      if(Flow % heat_transfer) then
        write(fu,'(a)') '#  1:y+  2:u+  3:kin_res+  4:kin_mod+  5:kin_tot+'  &
                     // '  6:uw_res+  7:uw_mod+  8:uw_tot+  9:vis_t  10:t+'  &
                     // '  11:t2+  12:ut+  13:vt+  14:wt+'
        do i = 1, n_prob
          if(n_count(i) .ne. 0) then
            write(fu,'(14es15.5e3)')                                  &
              wall_p(i),                                              &  !  1
              u_p(i),                                                 &  !  2
              0.5*(uu_p(i)+vv_p(i)+ww_p(i)),                          &  !  3
              kin_p(i),                                               &  !  4
              0.5*(uu_p(i)+vv_p(i)+ww_p(i)) + kin_p(i),               &  !  5
              uw_p(i),                                                &  !  6
              uw_mod_p(i),                                            &  !  7
              uw_p(i) + uw_mod_p(i),                                  &  !  8
              vis_t_p(i),                                             &  !  9
              t_p(i),                                                 &  ! 10
              t2_p(i),                                                &  ! 11
              ut_p(i),                                                &  ! 12
              vt_p(i),                                                &  ! 13
              wt_p(i)                                                    ! 14
          end if
        end do
      else
        write(fu,'(a)') '#  1:y+  2:u+  3:kin_res+  4:kin_mod+  5:kin_tot+'  &
                     // '  6:uw_res+  7:uw_mod+  8:uw_tot+  9:vis_t'
        do i = 1, n_prob
          if(n_count(i) .ne. 0) then
            write(fu,'(9es15.5e3)')                                   &
              wall_p(i),                                              &  !  1
              u_p(i),                                                 &  !  2
              0.5*(uu_p(i)+vv_p(i)+ww_p(i)),                          &  !  3
              kin_p(i),                                               &  !  4
              0.5*(uu_p(i)+vv_p(i)+ww_p(i)) + kin_p(i),               &  !  5
              uw_p(i),                                                &  !  6
              uw_mod_p(i),                                            &  !  7
              uw_p(i) + uw_mod_p(i),                                  &  !  8
              vis_t_p(i)                                                 !  9
          end if
        end do
      end if  ! heat_transfer

    else if(Turb % model == DES_SPALART) then

      !--- DES Spalart: no modelled k, but vis_t is available ------!
      if(Flow % heat_transfer) then
        write(fu,'(a)') '#  1:y+  2:u+  3:kin_res+  4:uw_res+  5:uw_mod+'    &
                     // '  6:uw_tot+  7:vis_t  8:t+  9:t2+  10:ut+  11:vt+'  &
                     // '  12:wt+'
        do i = 1, n_prob
          if(n_count(i) .ne. 0) then
            write(fu,'(12es15.5e3)')                                  &
              wall_p(i),                                              &  !  1
              u_p(i),                                                 &  !  2
              0.5*(uu_p(i)+vv_p(i)+ww_p(i)),                          &  !  3
              uw_p(i),                                                &  !  4
              uw_mod_p(i),                                            &  !  5
              uw_p(i) + uw_mod_p(i),                                  &  !  6
              vis_t_p(i),                                             &  !  7
              t_p(i),                                                 &  !  8
              t2_p(i),                                                &  !  9
              ut_p(i),                                                &  ! 10
              vt_p(i),                                                &  ! 11
              wt_p(i)                                                    ! 12
          end if
        end do
      else
        write(fu,'(a)') '#  1:y+  2:u+  3:kin_res+  4:uw_res+'       &
                     // '  5:uw_mod+  6:uw_tot+  7:vis_t'
        do i = 1, n_prob
          if(n_count(i) .ne. 0) then
            write(fu,'(7es15.5e3)')                                  &
              wall_p(i), u_p(i),                                     &  !  1-2
              0.5*(uu_p(i)+vv_p(i)+ww_p(i)),                         &  !  3
              uw_p(i), uw_mod_p(i), uw_p(i)+uw_mod_p(i),             &  !  4-6
              vis_t_p(i)                                                !  7
          end if
        end do
      end if

    else if(Turb % model == DNS              .or.  &
            Turb % model == NO_TURBULENCE_MODEL) then

      !--- DNS / no turbulence model: only resolved stresses --------!
      if(Flow % heat_transfer) then
        write(fu,'(a)') '#  1:y+  2:u+  3:kin_res+  4:uw_res+'              &
                     // '  5:t+  6:t2+  7:ut+  8:vt+  9:wt+'
        do i = 1, n_prob
          if(n_count(i) .ne. 0) then
            write(fu,'(9es15.5e3)')                                  &
              wall_p(i), u_p(i),                                     &  !  1-2
              0.5*(uu_p(i)+vv_p(i)+ww_p(i)),                         &  !  3
              uw_p(i),                                               &  !  4
              t_p(i), t2_p(i), ut_p(i), vt_p(i), wt_p(i)                !  5-9
          end if
        end do
      else
        write(fu,'(a)') '#  1:y+  2:u+  3:kin_res+  4:uw_res+'
        do i = 1, n_prob
          if(n_count(i) .ne. 0) then
            write(fu,'(4es15.5e3)')                                  &
              wall_p(i), u_p(i),                                     &  !  1-2
              0.5*(uu_p(i)+vv_p(i)+ww_p(i)),                         &  !  3
              uw_p(i)                                                   !  4
          end if
        end do
      end if

    else

      !--- Generic LES sub-grid model (Smagorinsky, Dynamic, WALE, !
      !   TVM, Prandtl-mixing): vis_t available, no modelled k     !
      if(Flow % heat_transfer) then
        write(fu,'(a)') '#  1:y+  2:u+  3:kin_res+  4:uw_res+'              &
                     // '  5:uw_mod+  6:uw_tot+  7:vis_t'                   &
                     // '  8:t+  9:t2+  10:ut+  11:vt+  12:wt+'
        do i = 1, n_prob
          if(n_count(i) .ne. 0) then
            write(fu,'(12es15.5e3)')                                 &
              wall_p(i), u_p(i),                                     &  !  1-2
              0.5*(uu_p(i)+vv_p(i)+ww_p(i)),                         &  !  3
              uw_p(i), uw_mod_p(i), uw_p(i)+uw_mod_p(i),             &  !  4-6
              vis_t_p(i),                                            &  !  7
              t_p(i), t2_p(i), ut_p(i), vt_p(i), wt_p(i)                !  8-12
          end if
        end do
      else
        write(fu,'(a)') '#  1:y+  2:u+  3:kin_res+  4:uw_res+'       &
                     // '  5:uw_mod+  6:uw_tot+  7:vis_t'
        do i = 1, n_prob
          if(n_count(i) .ne. 0) then
            write(fu,'(7es15.5e3)')                                  &
              wall_p(i), u_p(i),                                     &  !  1-2
              0.5*(uu_p(i)+vv_p(i)+ww_p(i)),                         &  !  3
              uw_p(i), uw_mod_p(i), uw_p(i)+uw_mod_p(i),             &  !  4-6
              vis_t_p(i)                                                !  7
          end if
        end do
      end if

    end if  ! model branch (use_stats)

  else

    !================================================================!
    !   RANS models (or LES/HYB in transient before stats window)   !
    !================================================================!
    select case(Turb % model)

      !--------------------------------------------------------------!
      case(K_EPS_ZETA_F)
      !--------------------------------------------------------------!
        if(Flow % heat_transfer) then
          write(fu,'(a)') '#  1:y+  2:u+  3:kin+  4:eps+  5:f22+'          &
                       // '  6:Re_t  7:vis_t  8:t+'
          do i = 1, n_prob
            if(n_count(i) .ne. 0) then
              write(fu,'(8es15.5e3)')                                      &
                wall_p(i), u_p(i), kin_p(i), eps_p(i),                     &
                f22_p(i), zeta_p(i), vis_t_p(i), t_p(i)
            end if
          end do
        else
          write(fu,'(a)') '#  1:y+  2:u+  3:kin+  4:eps+  5:f22+'          &
                       // '  6:Re_t  7:vis_t'
          do i = 1, n_prob
            if(n_count(i) .ne. 0) then
              write(fu,'(7es15.5e3)')                                      &
                wall_p(i), u_p(i), kin_p(i), eps_p(i),                     &
                f22_p(i), zeta_p(i), vis_t_p(i)
            end if
          end do
        end if

      !--------------------------------------------------------------!
      case(K_EPS)
      !--------------------------------------------------------------!
        if(Flow % heat_transfer) then
          write(fu,'(a)') '#  1:y+  2:u+  3:kin+  4:eps+  5:uw+  6:vis_t  7:t+'
          do i = 1, n_prob
            if(n_count(i) .ne. 0) then
              write(fu,'(7es15.5e3)')                                         &
                wall_p(i), u_p(i), kin_p(i), eps_p(i),                        &
                uw_p(i), vis_t_p(i), t_p(i)
            end if
          end do
        else
          write(fu,'(a)') '#  1:y+  2:u+  3:kin+  4:eps+  5:uw+  6:vis_t'
          do i = 1, n_prob
            if(n_count(i) .ne. 0) then
              write(fu,'(6es15.5e3)')                                         &
                wall_p(i), u_p(i), kin_p(i), eps_p(i),                        &
                uw_p(i), vis_t_p(i)
            end if
          end do
        end if

      !--------------------------------------------------------------!
      case(K_OMEGA_SST)
      !--------------------------------------------------------------!
        if(Flow % heat_transfer) then
          write(fu,'(a)') '#  1:y+  2:u+  3:kin+  4:omega+  5:uw+  6:vis_t  7:t+'
          do i = 1, n_prob
            if(n_count(i) .ne. 0) then
              write(fu,'(7es15.5e3)')                                         &
                wall_p(i), u_p(i), kin_p(i), omg_p(i),                        &
                uw_p(i), vis_t_p(i), t_p(i)
            end if
          end do
        else
          write(fu,'(a)') '#  1:y+  2:u+  3:kin+  4:omega+  5:uw+  6:vis_t'
          do i = 1, n_prob
            if(n_count(i) .ne. 0) then
              write(fu,'(6es15.5e3)')                                         &
                wall_p(i), u_p(i), kin_p(i), omg_p(i),                        &
                uw_p(i), vis_t_p(i)
            end if
          end do
        end if

      !--------------------------------------------------------------!
      case default  ! SPALART_ALLMARAS or any unrecognised RANS model
      !--------------------------------------------------------------!
        if(Flow % heat_transfer) then
          write(fu,'(a)') '#  1:y+  2:u+  3:uw+  4:vis_t  5:t+'
          do i = 1, n_prob
            if(n_count(i) .ne. 0) then
              write(fu,'(5es15.5e3)')                                         &
                wall_p(i), u_p(i), uw_p(i), vis_t_p(i), t_p(i)
            end if
          end do
        else
          write(fu,'(a)') '#  1:y+  2:u+  3:uw+  4:vis_t'
          do i = 1, n_prob
            if(n_count(i) .ne. 0) then
              write(fu,'(4es15.5e3)')                                         &
                wall_p(i), u_p(i), uw_p(i), vis_t_p(i)
            end if
          end do
        end if

    end select

  end if  ! use_stats / RANS output

  close(fu)

  end if  ! First_Proc

  if(First_Proc()) write(6,'(a)') '# Finished with User_Mod_Save_Results.'

  end subroutine
