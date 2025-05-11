!==============================================================================!
  subroutine Calculate_Mean(Turb, n0)
!------------------------------------------------------------------------------!
!   Calculates time averaged velocity and velocity fluctuations.               !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Turb_Type), target :: Turb
  integer,      intent(in) :: n0
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type), pointer :: Flow
  type(Grid_Type),  pointer :: Grid
  type(Var_Type),   pointer :: u, v, w, p, t, phi
  type(Var_Type),   pointer :: kin, eps, f22, zeta, vis, t2
  type(Var_Type),   pointer :: uu, vv, ww, uv, uw, vw
  type(Var_Type),   pointer :: ut, vt, wt
  integer                   :: c, n, sc
  real, contiguous, pointer :: u_mean(:), v_mean(:), w_mean(:),  &
                               p_mean(:), t_mean(:), q_mean(:)
  real, contiguous, pointer :: kin_mean (:), eps_mean(:),  &
                               zeta_mean(:), f22_mean(:)
  real, contiguous, pointer :: uu_res(:), vv_res(:), ww_res(:),  &
                               uv_res(:), vw_res(:), uw_res(:)
  real, contiguous, pointer :: ut_res(:), vt_res(:), wt_res(:), t2_res(:)
  real, contiguous, pointer :: ut_mean(:), vt_mean(:), wt_mean(:), t2_mean(:)
  real, contiguous, pointer :: phi_mean(:,:)
!==============================================================================!

  if(.not. Turb % statistics) return

  ! Take aliases
  Flow => Turb % pnt_flow
  Grid => Flow % pnt_grid
  p    => Flow % p
  vis  => Turb % vis
  t2   => Turb % t2
  call Flow % Alias_Momentum(u, v, w)
  call Flow % Alias_Energy  (t)
  call Turb % Alias_K_Eps_Zeta_F(kin, eps, zeta, f22)
  call Turb % Alias_Stresses    (uu, vv, ww, uv, uw, vw)
  call Turb % Alias_Heat_Fluxes (ut, vt, wt)

  ! Time averaged momentum and energy equations
  u_mean => Turb % u_mean;  v_mean => Turb % v_mean;  w_mean => Turb % w_mean
  p_mean => Turb % p_mean;  t_mean => Turb % t_mean;  q_mean => Turb % q_mean

  ! Time averaged modeled quantities
  if(Turb % model .eq. K_EPS) then
    kin_mean  => Turb % kin_mean;   eps_mean  => Turb % eps_mean
    if(Flow % heat_transfer) then
      ut_mean => Turb % ut_mean;  vt_mean => Turb % vt_mean
      wt_mean => Turb % wt_mean;  t2_mean => Turb % t2_mean
    end if
  end if

  ! Time averaged modeled quantities
  if(Turb % model .eq. K_EPS_ZETA_F  .or.  &
     Turb % model .eq. HYBRID_LES_RANS) then
    kin_mean  => Turb % kin_mean;   eps_mean  => Turb % eps_mean
    zeta_mean => Turb % zeta_mean;  f22_mean  => Turb % f22_mean
    if(Flow % heat_transfer) then
      ut_mean => Turb % ut_mean;  vt_mean => Turb % vt_mean
      wt_mean => Turb % wt_mean;  t2_mean => Turb % t2_mean
    end if
  end if

  ! Resolved Reynolds stresses and heat fluxes
  uu_res => Turb % uu_res;  vv_res => Turb % vv_res;  ww_res => Turb % ww_res
  uv_res => Turb % uv_res;  vw_res => Turb % vw_res;  uw_res => Turb % uw_res
  ut_res => Turb % ut_res;  vt_res => Turb % vt_res;  wt_res => Turb % wt_res
  t2_res => Turb % t2_res

  n = Time % Curr_Dt() - n0

  if(n > -1) then

    do c = Cells_At_Boundaries_In_Domain_And_Buffers()

      ! Mean velocities (and temperature)
      u_mean(c) = (u_mean(c) * real(n) + u % n(c)) / real(n+1)
      v_mean(c) = (v_mean(c) * real(n) + v % n(c)) / real(n+1)
      w_mean(c) = (w_mean(c) * real(n) + w % n(c)) / real(n+1)
      p_mean(c) = (p_mean(c) * real(n) + p % n(c)) / real(n+1)

      if(Flow % heat_transfer) then
        t_mean(c) = (t_mean(c) * real(n) + t % n(c)) / real(n+1)
        q_mean(c) = (q_mean(c) * real(n) + t % q(c)) / real(n+1)
      end if

      ! Resolved Reynolds stresses
      uu_res(c) = (uu_res(c) * real(n) + u % n(c) * u % n(c)) / real(n+1)
      vv_res(c) = (vv_res(c) * real(n) + v % n(c) * v % n(c)) / real(n+1)
      ww_res(c) = (ww_res(c) * real(n) + w % n(c) * w % n(c)) / real(n+1)

      uv_res(c) = (uv_res(c) * real(n) + u % n(c) * v % n(c)) / real(n+1)
      uw_res(c) = (uw_res(c) * real(n) + u % n(c) * w % n(c)) / real(n+1)
      vw_res(c) = (vw_res(c) * real(n) + v % n(c) * w % n(c)) / real(n+1)

      ! Resolved turbulent heat fluxes
      if(Flow % heat_transfer) then
        t2_res(c) = (t2_res(c) * real(n) + t % n(c) * t % n(c)) / real(n+1)
        ut_res(c) = (ut_res(c) * real(n) + u % n(c) * t % n(c)) / real(n+1)
        vt_res(c) = (vt_res(c) * real(n) + v % n(c) * t % n(c)) / real(n+1)
        wt_res(c) = (wt_res(c) * real(n) + w % n(c) * t % n(c)) / real(n+1)

        ! Resolved turbulent heat fluxes
        if(Turb % model .eq. K_EPS        .or.  &
           Turb % model .eq. K_EPS_ZETA_F .or.  &
           Turb % model .eq. HYBRID_LES_RANS) then

          t2_mean(c) = (t2_mean(c) * real(n) + t2 % n(c)) / real(n+1)
          ut_mean(c) = (ut_mean(c) * real(n) + ut % n(c)) / real(n+1)
          vt_mean(c) = (vt_mean(c) * real(n) + vt % n(c)) / real(n+1)
          wt_mean(c) = (wt_mean(c) * real(n) + wt % n(c)) / real(n+1)
        end if
      end if

      !-----------------!
      !   K-eps model   !
      !-----------------!
      if(Turb % model .eq. K_EPS) then

        ! Time-averaged modeled quantities
        kin_mean(c) = (kin_mean(c) * real(n) + kin % n(c)) / real(n+1)
        eps_mean(c) = (eps_mean(c) * real(n) + eps % n(c)) / real(n+1)
      end if

      !------------------!
      !   K-eps-zeta-f   !
      !------------------!
      if(Turb % model .eq. K_EPS_ZETA_F .or.  &
         Turb % model .eq. HYBRID_LES_RANS) then

        ! Time-averaged modeled quantities
        kin_mean (c) = (kin_mean (c) * real(n) + kin  % n(c)) / real(n+1)
        eps_mean (c) = (eps_mean (c) * real(n) + eps  % n(c)) / real(n+1)
        zeta_mean(c) = (zeta_mean(c) * real(n) + zeta % n(c)) / real(n+1)
        f22_mean (c) = (f22_mean (c) * real(n) + f22  % n(c)) / real(n+1)
      end if

      !-------------!
      !   Scalars   !
      !-------------!
      do sc = 1, Flow % n_scalars
        phi      => Flow % scalar(sc)
        phi_mean => Turb % scalar_mean
        phi_mean(sc, c) = (phi_mean(sc, c) * real(n) + phi % n(c)) / real(n+1)
      end do
    end do

  end if

  end subroutine
