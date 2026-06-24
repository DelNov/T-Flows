!==============================================================================!
  subroutine User_Mod_Save_Results(Flow, Turb, Vof, Swarm, domain)
!------------------------------------------------------------------------------!
!   Common channel-flow post-processing routine for RANS, LES and HYB models.   !
!                                                                              !
!   For RANS models, instantaneous fields are used. For LES/HYB models, mean    !
!   fields and resolved stresses are used only after statistics have started    !
!   and the corresponding arrays are allocated. Before that, instantaneous      !
!   fields are used to avoid accessing unallocated or empty statistics arrays.  !
!                                                                              !
!   The Dean estimate of Cf is replaced by the Zanoun et al. channel-flow       !
!   correlation. Dean's formula is used only as the initial Newton guess.       !
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
  type(Var_Type),  pointer :: u, v, w, kin, eps, zeta, f22, t
  integer                  :: n_prob, pl, c, i, count, n_points
  integer                  :: fu, c1, c2, s, iter
  character(SL)            :: res_name
  real, allocatable        :: wall_p(:), u_p(:), v_p(:), w_p(:),        &
                              y_plus_p(:), kin_mod_p(:), kin_res_p(:), &
                              uw_mod_p(:), uw_res_p(:), vis_t_p(:),    &
                              t_p(:)
  integer, allocatable     :: n_count(:)
  real                     :: ubulk, error, re, cf_zanoun, cf
  real                     :: u_tau_p, t_tau, t_wall
  real                     :: dens_const, visc_const
  real                     :: capa_const, cond_const
  real                     :: u_loc, v_loc, w_loc, t_loc
  real                     :: kin_mod_loc, kin_res_loc, uw_res_loc
  real                     :: uw_mod_loc, vis_t_loc
  real                     :: zan_fun, zan_der
  logical                  :: use_mean, use_mean_heat
  logical                  :: les_hyb_model
  logical                  :: k_eps_like, k_omega_like, has_vis_t
!==============================================================================!

  ! Don't save if this is initial condition, nothing is developed yet
  if(Time % Curr_Dt() .eq. 0) return

  ! Take aliases
  Grid => Flow % pnt_grid
  bulk => Flow % bulk
  call Flow % Alias_Momentum(u, v, w)
  call Flow % Alias_Energy  (t)
  nullify(kin, eps, zeta, f22)

  ! Take constant physical properties
  call Control % Mass_Density        (dens_const)
  call Control % Dynamic_Viscosity   (visc_const)
  call Control % Heat_Capacity       (capa_const)
  call Control % Thermal_Conductivity(cond_const)

  ! Set file name to save results
  call File % Set_Name(res_name,                      &
                       time_step = Time % Curr_Dt(),  &
                       extension = '-res.dat')

  ! Number of probes in cell rows
  n_prob = Grid % n_z_planes - 1

  ubulk    = bulk % flux_x / bulk % area_x
  t_wall   = 0.0
  n_points = 0

  ! LES/HYB statistics should be used only after they have started.
  les_hyb_model = Turb % model == LES_SMAGORINSKY  .or.  &
                  Turb % model == LES_WALE         .or.  &
                  Turb % model == LES_DYNAMIC      .or.  &
                  Turb % model == LES_TVM          .or.  &
                  Turb % model == HYBRID_LES_PRANDTL .or. &
                  Turb % model == HYBRID_LES_RANS  .or.  &
                  Turb % model == DES_SPALART

  use_mean = .false.
  if(les_hyb_model) then
    if(Turb % statistics) then
      if(allocated(Turb % u_mean) .and.  &
         allocated(Turb % v_mean) .and.  &
         allocated(Turb % w_mean) .and.  &
         allocated(Turb % uu_res) .and.  &
         allocated(Turb % vv_res) .and.  &
         allocated(Turb % ww_res) .and.  &
         allocated(Turb % uw_res)) then
        use_mean = .true.
      end if
    end if
  end if

  use_mean_heat = .false.
  if(use_mean .and. Flow % heat_transfer) then
    if(allocated(Turb % t_mean)) use_mean_heat = .true.
  end if

  k_eps_like = Turb % model == K_EPS         .or.  &
               Turb % model == K_EPS_ZETA_F .or.  &
               Turb % model == HYBRID_LES_RANS

  k_omega_like = Turb % model == K_OMEGA_SST

  if(k_eps_like .or. k_omega_like) then
    call Turb % Alias_K_Eps_Zeta_F(kin, eps, zeta, f22)
  end if

  has_vis_t = Turb % model == K_EPS            .or.  &
              Turb % model == K_EPS_ZETA_F    .or.  &
              Turb % model == K_OMEGA_SST     .or.  &
              Turb % model == SPALART_ALLMARAS.or.  &
              Turb % model == DES_SPALART     .or.  &
              Turb % model == LES_SMAGORINSKY .or.  &
              Turb % model == LES_WALE        .or.  &
              Turb % model == LES_DYNAMIC     .or.  &
              Turb % model == LES_TVM         .or.  &
              Turb % model == HYBRID_LES_PRANDTL .or. &
              Turb % model == HYBRID_LES_RANS

  !-------------------------------------------------------------------------!
  !   Allocate memory for variables to be extracted in homogeneous planes   !
  !-------------------------------------------------------------------------!
  allocate(wall_p    (n_prob));  wall_p     = 0.0
  allocate(u_p       (n_prob));  u_p        = 0.0
  allocate(v_p       (n_prob));  v_p        = 0.0
  allocate(w_p       (n_prob));  w_p        = 0.0
  allocate(kin_mod_p (n_prob));  kin_mod_p  = 0.0
  allocate(kin_res_p (n_prob));  kin_res_p  = 0.0
  allocate(uw_mod_p  (n_prob));  uw_mod_p   = 0.0
  allocate(uw_res_p  (n_prob));  uw_res_p   = 0.0
  allocate(vis_t_p   (n_prob));  vis_t_p    = 0.0
  allocate(y_plus_p  (n_prob));  y_plus_p   = 0.0
  allocate(n_count   (n_prob));  n_count    = 0

  if(Flow % heat_transfer) then
    allocate(t_p(n_prob));  t_p = 0.0
  end if

  count = 0

  !-------------------------!
  !   Average the results   !
  !-------------------------!
  do i = 1, n_prob
    do c = Cells_In_Domain()
      if(Grid % zc(c) > Grid % z_coord_plane(i) .and.  &
         Grid % zc(c) < Grid % z_coord_plane(i+1)) then

        if(use_mean) then
          u_loc = Turb % u_mean(c)
          v_loc = Turb % v_mean(c)
          w_loc = Turb % w_mean(c)

          kin_res_loc = 0.5 * (Turb % uu_res(c)  &
                              - Turb % u_mean(c) * Turb % u_mean(c)  &
                              + Turb % vv_res(c)  &
                              - Turb % v_mean(c) * Turb % v_mean(c)  &
                              + Turb % ww_res(c)  &
                              - Turb % w_mean(c) * Turb % w_mean(c))

          uw_res_loc = Turb % uw_res(c)  &
                     - Turb % u_mean(c) * Turb % w_mean(c)
        else
          u_loc       = u % n(c)
          v_loc       = v % n(c)
          w_loc       = w % n(c)
          kin_res_loc = 0.0
          uw_res_loc  = 0.0
        end if

        kin_mod_loc = 0.0
        if(use_mean .and. allocated(Turb % kin_mean)) then
          kin_mod_loc = Turb % kin_mean(c)
        else if(k_eps_like) then
          kin_mod_loc = kin % n(c)
        else if(k_omega_like) then
          kin_mod_loc = kin % n(c)
        end if

        vis_t_loc = 0.0
        if(Turb % model == HYBRID_LES_RANS) then
          vis_t_loc = Turb % vis_t_eff(c)
        else if(has_vis_t) then
          vis_t_loc = Turb % vis_t(c)
        end if

        uw_mod_loc = vis_t_loc * (u % z(c) + w % x(c))

        wall_p    (i) = wall_p    (i) + Grid % wall_dist(c)
        u_p       (i) = u_p       (i) + u_loc
        v_p       (i) = v_p       (i) + v_loc
        w_p       (i) = w_p       (i) + w_loc
        kin_res_p (i) = kin_res_p (i) + kin_res_loc
        kin_mod_p (i) = kin_mod_p (i) + kin_mod_loc
        uw_res_p  (i) = uw_res_p  (i) + uw_res_loc
        uw_mod_p  (i) = uw_mod_p  (i) + uw_mod_loc
        vis_t_p   (i) = vis_t_p   (i) + vis_t_loc / visc_const
        y_plus_p  (i) = y_plus_p  (i) + Turb % y_plus(c)

        if(Flow % heat_transfer) then
          if(use_mean_heat) then
            t_loc = Turb % t_mean(c)
          else
            t_loc = t % n(c)
          end if
          t_p(i) = t_p(i) + t_loc
        end if

        n_count(i) = n_count(i) + 1
      end if
    end do
  end do

  ! Average over all processors
  do pl = 1, n_prob
    call Global % Sum_Int (n_count(pl))
    call Global % Sum_Real(wall_p   (pl))
    call Global % Sum_Real(u_p      (pl))
    call Global % Sum_Real(v_p      (pl))
    call Global % Sum_Real(w_p      (pl))
    call Global % Sum_Real(kin_res_p(pl))
    call Global % Sum_Real(kin_mod_p(pl))
    call Global % Sum_Real(uw_res_p (pl))
    call Global % Sum_Real(uw_mod_p (pl))
    call Global % Sum_Real(vis_t_p  (pl))
    call Global % Sum_Real(y_plus_p (pl))

    count = count + n_count(pl)

    if(Flow % heat_transfer) then
      call Global % Sum_Real(t_p(pl))
    end if
  end do

  call Global % Wait

  do i = 1, n_prob
    if(n_count(i) .ne. 0) then
      wall_p    (i) = wall_p    (i) / n_count(i)
      u_p       (i) = u_p       (i) / n_count(i)
      v_p       (i) = v_p       (i) / n_count(i)
      w_p       (i) = w_p       (i) / n_count(i)
      kin_res_p (i) = kin_res_p (i) / n_count(i)
      kin_mod_p (i) = kin_mod_p (i) / n_count(i)
      uw_res_p  (i) = uw_res_p  (i) / n_count(i)
      uw_mod_p  (i) = uw_mod_p  (i) / n_count(i)
      vis_t_p   (i) = vis_t_p   (i) / n_count(i)
      y_plus_p  (i) = y_plus_p  (i) / n_count(i)

      if(Flow % heat_transfer) then
        t_p(i) = t_p(i) / n_count(i)
      end if
    end if
  end do

  !------------------------------------!
  !   Non-dimensionalize the results   !
  !------------------------------------!

  ! Calculating friction velocity
  if(y_plus_p(1) > 5.0) then
    u_tau_p = sqrt(max(abs(bulk % p_drop_x),  &
                       abs(bulk % p_drop_y),  &
                       abs(bulk % p_drop_z)) / dens_const)
  else
    u_tau_p = sqrt((visc_const * sqrt(u_p(1)**2 + v_p(1)**2 + w_p(1)**2)  &
                  / wall_p(1)) / dens_const)
  end if

  if(u_tau_p .eq. 0.0) then
    if(First_Proc()) then
      write(*,*) '# Friction velocity is zero in User_Mod_Save_Results.f90!'
    end if
    return
  end if

  if(Flow % heat_transfer) then
    t_tau = Flow % heat_flux / (dens_const * capa_const * u_tau_p)

    do s = 1, Grid % n_faces
      c1 = Grid % faces_c(1,s)
      c2 = Grid % faces_c(2,s)
      if(c2 < 0) then
        if(Grid % Bnd_Cond_Type(c2) .eq. WALL .or.  &
           Grid % Bnd_Cond_Type(c2) .eq. WALLFL) then

          if(use_mean_heat) then
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

  !----------------------------------------!
  !   Write the results in the .dat file   !
  !----------------------------------------!
  call File % Open_For_Writing_Ascii(res_name, fu)

  re = abs(dens_const * ubulk * 2.0 / visc_const)

  ! Zanoun et al. skin-friction correlation for plane channel flow:
  !   1/sqrt(Cf) = 1.911*ln(Re*sqrt(Cf)) - 1.282
  ! The equation is implicit in Cf and is solved here by Newton iteration.
  cf_zanoun = 0.073 * re**(-0.25)   ! Dean correlation used only as initial guess

  do iter = 1, 50
    zan_fun = 1.0 / sqrt(cf_zanoun)                         &
            - 1.911 * log(re * sqrt(cf_zanoun)) + 1.282

    zan_der = -0.5 * cf_zanoun**(-1.5)                      &
              - 0.5 * 1.911 / cf_zanoun

    cf_zanoun = max(1.0e-12, cf_zanoun - zan_fun / zan_der)

    if(abs(zan_fun) < 1.0e-8) exit
  end do

  cf    = u_tau_p**2 / (0.5 * ubulk**2)
  error = abs(cf_zanoun - cf) / cf_zanoun * 100.0

  write(fu,'(a1,(a12,e12.6))')  &
  '#', 'Density  = ', dens_const
  write(fu,'(a1,(a12,e12.6))')  &
  '#', 'Ubulk    = ', ubulk
  write(fu,'(a1,(a12,e12.6))')  &
  '#', 'Re       = ', re
  write(fu,'(a1,(a12,e12.6))')  &
  '#', 'Re_tau   = ', dens_const * u_tau_p / visc_const
  write(fu,'(a1,(a12,e12.6))')  &
  '#', 'Cf       = ', 2.0 * (u_tau_p / ubulk)**2
  write(fu,'(a1,(a12,f12.6))')  &
  '#', 'Utau     = ', u_tau_p
  write(fu,'(a1,(a12,f12.6,a2,a22))') &
  '#', 'Cf_error = ', error, ' %', 'Zanoun formula used.'
  write(fu,'(a1,(a12,l1))') &
  '#', 'Mean used= ', use_mean

  do i = 1, n_prob
    wall_p    (i) = dens_const * wall_p(i) * u_tau_p / visc_const
    u_p       (i) = u_p      (i) / u_tau_p
    v_p       (i) = v_p      (i) / u_tau_p
    w_p       (i) = w_p      (i) / u_tau_p
    kin_res_p (i) = kin_res_p(i) / u_tau_p**2
    kin_mod_p (i) = kin_mod_p(i) / u_tau_p**2
    uw_res_p  (i) = uw_res_p (i) / (u_tau_p**2 * dens_const)
    uw_mod_p  (i) = uw_mod_p (i) / (u_tau_p**2 * dens_const)

    if(Flow % heat_transfer) then
      if(abs(t_tau) > TINY) then
        t_p(i) = (t_wall - t_p(i)) / t_tau
      else
        t_p(i) = 0.0
      end if
    end if
  end do

  if(Flow % heat_transfer) then
    write(fu,'(a)') '#  1:z+,  2:u+,  3:k_res+,  4:k_mod+,'    //  &
                    '  5:k_tot+,  6:uw_res+,  7:uw_mod+,'      //  &
                    '  8:uw_tot+,  9:vis_t/visc_const,  10:T+'
    do i = 1, n_prob
      if(n_count(i) .ne. 0) then
        write(fu,'(10es15.5e3)') wall_p(i),                    &
                                  u_p(i),                       &
                                  kin_res_p(i),                 &
                                  kin_mod_p(i),                 &
                                  kin_res_p(i) + kin_mod_p(i),  &
                                  uw_res_p(i),                  &
                                  uw_mod_p(i),                  &
                                  uw_res_p(i) + uw_mod_p(i),    &
                                  vis_t_p(i),                   &
                                  t_p(i)
      end if
    end do
  else
    write(fu,'(a)') '#  1:z+,  2:u+,  3:k_res+,  4:k_mod+,'    //  &
                    '  5:k_tot+,  6:uw_res+,  7:uw_mod+,'      //  &
                    '  8:uw_tot+,  9:vis_t/visc_const'
    do i = 1, n_prob
      if(n_count(i) .ne. 0) then
        write(fu,'(9es15.5e3)')  wall_p(i),                    &
                                  u_p(i),                       &
                                  kin_res_p(i),                 &
                                  kin_mod_p(i),                 &
                                  kin_res_p(i) + kin_mod_p(i),  &
                                  uw_res_p(i),                  &
                                  uw_mod_p(i),                  &
                                  uw_res_p(i) + uw_mod_p(i),    &
                                  vis_t_p(i)
      end if
    end do
  end if

  close(fu)

  deallocate(wall_p)
  deallocate(u_p)
  deallocate(v_p)
  deallocate(w_p)
  deallocate(kin_res_p)
  deallocate(kin_mod_p)
  deallocate(uw_res_p)
  deallocate(uw_mod_p)
  deallocate(vis_t_p)
  deallocate(y_plus_p)
  deallocate(n_count)
  if(Flow % heat_transfer) deallocate(t_p)

  if(First_Proc()) write(6, *) '# Finished with User_Mod_Save_Results.f90.'

  end subroutine
