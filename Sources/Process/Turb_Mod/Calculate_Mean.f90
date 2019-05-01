!==============================================================================!
  subroutine Turb_Mod_Calculate_Mean(turb, n0, n1)
!------------------------------------------------------------------------------!
!   Calculates time averaged velocity and velocity fluctuations.               !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Turb_Type),  target :: turb
  integer                  :: n0, n1
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type), pointer :: flow
  type(Grid_Type),  pointer :: grid
  type(Var_Type),   pointer :: u, v, w, p, t
  type(Var_Type),   pointer :: kin, eps, f22, zeta, vis, t2
  type(Var_Type),   pointer :: uu, vv, ww, uv, uw, vw
  integer                   :: c, n
  real,             pointer :: u_mean(:), v_mean(:), w_mean(:),  &
                               p_mean(:), t_mean(:)
  real,             pointer :: kin_mean (:), eps_mean(:),  &
                               zeta_mean(:), f22_mean(:)
  real,             pointer :: uu_res(:), vv_res(:), ww_res(:),  &
                               uv_res(:), vw_res(:), uw_res(:)
  real,             pointer :: ut_res(:), vt_res(:), wt_res(:), t2_res(:)
  real,             pointer :: uu_mean(:), vv_mean(:), ww_mean(:)
  real,             pointer :: uv_mean(:), vw_mean(:), uw_mean(:)
  real,             pointer :: ut_mean(:), vt_mean(:), wt_mean(:), t2_mean(:)
!==============================================================================!

  if(.not. turbulence_statistics) return

  ! Take aliases
  flow => turb % pnt_flow
  grid => flow % pnt_grid
  p    => flow % p
  vis  => turb % vis
  t2   => turb % t2
  call Field_Mod_Alias_Momentum   (flow, u, v, w)
  call Field_Mod_Alias_Energy     (flow, t)
  call Turb_Mod_Alias_K_Eps_Zeta_F(turb, kin, eps, zeta, f22)
  call Turb_Mod_Alias_Stresses    (turb, uu, vv, ww, uv, uw, vw)

  ! Time averaged momentum and energy equations
  u_mean => turb % u_mean;  v_mean => turb % v_mean;  w_mean => turb % w_mean
  p_mean => turb % p_mean;  t_mean => turb % t_mean

  ! Time averaged modeled quantities
  kin_mean  => turb % kin_mean;   eps_mean  => turb % eps_mean
  zeta_mean => turb % zeta_mean;  f22_mean  => turb % f22_mean

  ! Time-averaged modelled Reynolds stresses and heat fluxes
  uu_mean => turb % uu_mean;  vv_mean => turb % vv_mean
  ww_mean => turb % ww_mean;  uv_mean => turb % uv_mean
  vw_mean => turb % vw_mean;  uw_mean => turb % uw_mean
  ut_mean => turb % ut_mean;  vt_mean => turb % vt_mean
  wt_mean => turb % wt_mean;  t2_mean => turb % t2_mean

  ! Resolved Reynolds stresses and heat fluxes
  uu_res => turb % uu_res;  vv_res => turb % vv_res;  ww_res => turb % ww_res
  uv_res => turb % uv_res;  vw_res => turb % vw_res;  uw_res => turb % uw_res
  ut_res => turb % ut_res;  vt_res => turb % vt_res;  wt_res => turb % wt_res
  t2_res => turb % t2_res

  n = n1 - n0

  if(n > -1) then

    do c = -grid % n_bnd_cells, grid % n_cells

      !---------------------------------!
      !   Scale-resolving simulations   ! 
      !---------------------------------!
      if(turbulence_model .eq. LES_SMAGORINSKY .or.  &
         turbulence_model .eq. LES_DYNAMIC     .or.  &
         turbulence_model .eq. LES_WALE        .or.  &
         turbulence_model .eq. DES_SPALART     .or.  &
         turbulence_model .eq. DNS) then

        ! Mean velocities (and temperature)
        u_mean(c) = (u_mean(c) * (1.*n) + u % n(c)) / (1.*(n+1))
        v_mean(c) = (v_mean(c) * (1.*n) + v % n(c)) / (1.*(n+1))
        w_mean(c) = (w_mean(c) * (1.*n) + w % n(c)) / (1.*(n+1))
        p_mean(c) = (p_mean(c) * (1.*n) + p % n(c)) / (1.*(n+1))
        if(heat_transfer) then
          t_mean(c) = (t_mean(c) * (1.*n) + t % n(c)) / (1.*(n+1))
        end if

        ! Resolved Reynolds stresses
        uu_res(c) = (uu_res(c)*(1.*n) + u % n(c) * u % n(c)) / (1.*(n+1))
        vv_res(c) = (vv_res(c)*(1.*n) + v % n(c) * v % n(c)) / (1.*(n+1))
        ww_res(c) = (ww_res(c)*(1.*n) + w % n(c) * w % n(c)) / (1.*(n+1))

        uv_res(c) = (uv_res(c)*(1.*n) + u % n(c) * v % n(c)) / (1.*(n+1))
        uw_res(c) = (uw_res(c)*(1.*n) + u % n(c) * w % n(c)) / (1.*(n+1))
        vw_res(c) = (vw_res(c)*(1.*n) + v % n(c) * w % n(c)) / (1.*(n+1))

        ! Resolved turbulent heat fluxes
        if(heat_transfer) then
          t2_res(c) = (t2_res(c)*(1.*n) + t % n(c) * t % n(c)) / (1.*(n+1))
          ut_res(c) = (ut_res(c)*(1.*n) + u % n(c) * t % n(c)) / (1.*(n+1))
          vt_res(c) = (vt_res(c)*(1.*n) + v % n(c) * t % n(c)) / (1.*(n+1))
          wt_res(c) = (wt_res(c)*(1.*n) + w % n(c) * t % n(c)) / (1.*(n+1))
        end if
      end if

      !-----------------!
      !   K-eps model   !
      !-----------------!
      if(turbulence_model .eq. K_EPS) then

        ! Time-averaged velocities (and temperature)
        u_mean(c) = (u_mean(c) * (1.*n) + u % n(c)) / (1.*(n+1))
        v_mean(c) = (v_mean(c) * (1.*n) + v % n(c)) / (1.*(n+1))
        w_mean(c) = (w_mean(c) * (1.*n) + w % n(c)) / (1.*(n+1))
        p_mean(c) = (p_mean(c) * (1.*n) + p % n(c)) / (1.*(n+1))
        if(heat_transfer) then
          t_mean(c) = (t_mean(c) * (1.*n) + t % n(c)) / (1.*(n+1))
        end if

        ! Time-averaged modeled quantities
        kin_mean(c) = (kin_mean(c) * (1.*n) + kin % n(c)) / (1.*(n+1))
        eps_mean(c) = (eps_mean(c) * (1.*n) + eps % n(c)) / (1.*(n+1))

        ! Resolved Reynolds stresses
        uu_res(c) = (uu_res(c)*(1.*n) + u % n(c) * u % n(c)) / (1.*(n+1))
        vv_res(c) = (vv_res(c)*(1.*n) + v % n(c) * v % n(c)) / (1.*(n+1))
        ww_res(c) = (ww_res(c)*(1.*n) + w % n(c) * w % n(c)) / (1.*(n+1))

        uv_res(c) = (uv_res(c)*(1.*n) + u % n(c) * v % n(c)) / (1.*(n+1))
        uw_res(c) = (uw_res(c)*(1.*n) + u % n(c) * w % n(c)) / (1.*(n+1))
        vw_res(c) = (vw_res(c)*(1.*n) + v % n(c) * w % n(c)) / (1.*(n+1))

        ! Resolved turbulent heat fluxes
        if(heat_transfer) then
          t2_res(c) = (t2_res(c)*(1.*n) + t % n(c) * t % n(c)) / (1.*(n+1))
          ut_res(c) = (ut_res(c)*(1.*n) + u % n(c) * t % n(c)) / (1.*(n+1))
          vt_res(c) = (vt_res(c)*(1.*n) + v % n(c) * t % n(c)) / (1.*(n+1))
          wt_res(c) = (wt_res(c)*(1.*n) + w % n(c) * t % n(c)) / (1.*(n+1))
        end if
      end if

      !------------------!
      !   K-eps-zeta-f   !
      !------------------!
      if(turbulence_model .eq. K_EPS_ZETA_F .or.  &
         turbulence_model .eq. HYBRID_LES_RANS) then

        ! Time-averaged velocities (and temperature)
        u_mean(c) = (u_mean(c) * (1.*n) + u % n(c)) / (1.*(n+1))
        v_mean(c) = (v_mean(c) * (1.*n) + v % n(c)) / (1.*(n+1))
        w_mean(c) = (w_mean(c) * (1.*n) + w % n(c)) / (1.*(n+1))
        p_mean(c) = (p_mean(c) * (1.*n) + p % n(c)) / (1.*(n+1))
        if (heat_transfer) then
          t_mean(c) = (t_mean(c) * (1.*n) + t % n(c)) / (1.*(n+1))
        end if

        ! Time-averaged modeled quantities
        kin_mean (c) = (kin_mean (c) * (1.*n) + kin  % n(c)) / (1.*(n+1))
        eps_mean (c) = (eps_mean (c) * (1.*n) + eps  % n(c)) / (1.*(n+1))
        zeta_mean(c) = (zeta_mean(c) * (1.*n) + zeta % n(c)) / (1.*(n+1))
        f22_mean (c) = (f22_mean (c) * (1.*n) + f22  % n(c)) / (1.*(n+1))
        if (heat_transfer) then
          t2_mean(c) = (t2_mean(c) * (1.*n) + t2 % n(c)) / (1.*(n+1))
        end if

        ! Resolved Reynolds stresses
        uu_res(c) = (uu_res(c)*(1.*n) + u % n(c) * u % n(c)) / (1.*(n+1))
        vv_res(c) = (vv_res(c)*(1.*n) + v % n(c) * v % n(c)) / (1.*(n+1))
        ww_res(c) = (ww_res(c)*(1.*n) + w % n(c) * w % n(c)) / (1.*(n+1))

        uv_res(c) = (uv_res(c)*(1.*n) + u % n(c) * v % n(c)) / (1.*(n+1))
        uw_res(c) = (uw_res(c)*(1.*n) + u % n(c) * w % n(c)) / (1.*(n+1))
        vw_res(c) = (vw_res(c)*(1.*n) + v % n(c) * w % n(c)) / (1.*(n+1))

        ! Resolved turbulent heat fluxes
        if (heat_transfer) then
          t2_res(c) = (t2_res(c)*(1.*n) + t % n(c) * t % n(c)) / (1.*(n+1))
          ut_res(c) = (ut_res(c)*(1.*n) + u % n(c) * t % n(c)) / (1.*(n+1))
          vt_res(c) = (vt_res(c)*(1.*n) + v % n(c) * t % n(c)) / (1.*(n+1))
          wt_res(c) = (wt_res(c)*(1.*n) + w % n(c) * t % n(c)) / (1.*(n+1))
        end if

      end if

      !----------------------------!
      !   Reynolds stress models   !
      !----------------------------!
      if(turbulence_model .eq. RSM_HANJALIC_JAKIRLIC .or.  &
         turbulence_model .eq. RSM_MANCEAU_HANJALIC) then

        ! Time-averaged velocities (and temperature)
        u_mean(c) = (u_mean(c) * (1.*n) + u % n(c)) / (1.*(n+1))
        v_mean(c) = (v_mean(c) * (1.*n) + v % n(c)) / (1.*(n+1))
        w_mean(c) = (w_mean(c) * (1.*n) + w % n(c)) / (1.*(n+1))
        p_mean(c) = (p_mean(c) * (1.*n) + p % n(c)) / (1.*(n+1))
        if (heat_transfer) then
          t_mean(c) = (t_mean(c) * (1.*n) + t % n(c)) / (1.*(n+1))
        end if

        ! Time-averaged modeled quantities (modelled Reynolds stresses)
        uu_mean (c) = (uu_mean (c) * (1.*n) + uu  % n(c)) / (1.*(n+1))
        vv_mean (c) = (vv_mean (c) * (1.*n) + vv  % n(c)) / (1.*(n+1))
        ww_mean (c) = (ww_mean (c) * (1.*n) + ww  % n(c)) / (1.*(n+1))
        uv_mean (c) = (uv_mean (c) * (1.*n) + uv  % n(c)) / (1.*(n+1))
        uw_mean (c) = (uw_mean (c) * (1.*n) + uw  % n(c)) / (1.*(n+1))
        vw_mean (c) = (vw_mean (c) * (1.*n) + vw  % n(c)) / (1.*(n+1))
        kin_mean(c) = (kin_mean(c) * (1.*n) + kin % n(c)) / (1.*(n+1))
        eps_mean(c) = (eps_mean(c) * (1.*n) + eps % n(c)) / (1.*(n+1))
        if(turbulence_model .eq. RSM_MANCEAU_HANJALIC) then
          f22_mean(c) = (f22_mean(c) * (1.*n) + f22 % n(c)) / (1.*(n+1))
        end if

        ! Resolved Reynolds stresses
        uu_res(c) = (uu_res(c)*(1.*n) + u % n(c) * u % n(c)) / (1.*(n+1))
        vv_res(c) = (vv_res(c)*(1.*n) + v % n(c) * v % n(c)) / (1.*(n+1))
        ww_res(c) = (ww_res(c)*(1.*n) + w % n(c) * w % n(c)) / (1.*(n+1))

        uv_res(c) = (uv_res(c)*(1.*n) + u % n(c) * v % n(c)) / (1.*(n+1))
        uw_res(c) = (uw_res(c)*(1.*n) + u % n(c) * w % n(c)) / (1.*(n+1))
        vw_res(c) = (vw_res(c)*(1.*n) + v % n(c) * w % n(c)) / (1.*(n+1))

        ! Resolved turbulent heat fluxes
        if (heat_transfer) then
          t2_res(c) = (t2_res(c)*(1.*n) + t % n(c) * t % n(c)) / (1.*(n+1))
          ut_res(c) = (ut_res(c)*(1.*n) + u % n(c) * t % n(c)) / (1.*(n+1))
          vt_res(c) = (vt_res(c)*(1.*n) + v % n(c) * t % n(c)) / (1.*(n+1))
          wt_res(c) = (wt_res(c)*(1.*n) + w % n(c) * t % n(c)) / (1.*(n+1))
        end if

      end if

      !------------------------------!
      !   User scalars are missing   !
      !------------------------------!
    end do 
  end if

  end subroutine
