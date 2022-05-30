!==============================================================================!
  subroutine Calculate_Heat_Flux(Turb)
!------------------------------------------------------------------------------!
!   Computes turbulent heat fluxes                                             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Turb_Type), target :: Turb
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type), pointer :: Flow
  type(Grid_Type),  pointer :: Grid
  type(Var_Type),   pointer :: uu, vv, ww, uv, uw, vw
  type(Var_Type),   pointer :: t, ut, vt, wt, t2
  type(Var_Type),   pointer :: u, v, w
  integer                   :: c, k, c1, c2, s
  real                      :: ut_log_law, vt_log_law, wt_log_law
  real                      :: nx, ny, nz, qx, qy, qz, ebf
!==============================================================================!

  ! Take aliases
  Flow => Turb % pnt_flow
  Grid => Flow % pnt_grid
  t    => Flow % t
  call Turb % Alias_Heat_Fluxes(ut, vt, wt)
  call Turb % Alias_Stresses   (uu, vv, ww, uv, uw, vw)
  call Flow % Alias_Momentum   (u, v, w)
  call Turb % Alias_T2         (t2)

  ! Check if these are already computed somewhere, ...
  ! ... maybe this call is not needed
  call Flow % Grad_Variable(t)

  !-----------------------------------------!
  !   Compute the sources in the interior   !
  !-----------------------------------------!
  call Control_Mod_Turbulent_Prandtl_Number(pr_t)

  !-----------------------------------------------------------!
  !   By default turbulent heat flux is caalculated by SGDH   !
  !-----------------------------------------------------------!

  if(Turb % heat_flux_model .eq. SGDH) then
    do c = 1, Grid % n_cells
      pr_t = max(Turb % Prandtl_Turb(c), TINY)
      ut % n(c) = - Turb % vis_t(c) / Flow % density(c) / pr_t * t % x(c)
      vt % n(c) = - Turb % vis_t(c) / Flow % density(c) / pr_t * t % y(c)
      wt % n(c) = - Turb % vis_t(c) / Flow % density(c) / pr_t * t % z(c)

      if(Turb % model .eq. HYBRID_LES_RANS) then
        ut % n(c) = - Turb % vis_t_eff(c) / Flow % density(c) &
                                          / pr_t * t % x(c)
        vt % n(c) = - Turb % vis_t_eff(c) / Flow % density(c) &
                                          / pr_t * t % y(c)
        wt % n(c) = - Turb % vis_t_eff(c) / Flow % density(c) &
                                          / pr_t * t % z(c)
      end if
    end do

  else if(Turb % heat_flux_model .eq. GGDH) then
    do c = 1, Grid % n_cells
      ut % n(c) = -c_theta * Turb % t_scale(c) * (uu % n(c) * t % x(c)  +  &
                                                  uv % n(c) * t % y(c)  +  &
                                                  uw % n(c) * t % z(c))
      vt % n(c) = -c_theta * Turb % t_scale(c) * (uv % n(c) * t % x(c)  +  &
                                                  vv % n(c) * t % y(c)  +  &
                                                  vw % n(c) * t % z(c))
      wt % n(c) = -c_theta * Turb % t_scale(c) * (uw % n(c) * t % x(c)  +  &
                                                  vw % n(c) * t % y(c)  +  &
                                                  ww % n(c) * t % z(c))
    end do

  else if(Turb % heat_flux_model .eq. AFM) then
    call Flow % Grad_Variable(Flow % u)
    call Flow % Grad_Variable(Flow % v)
    call Flow % Grad_Variable(Flow % w)

    do c = 1, Grid % n_cells
      ut % n(c) = -c_theta * Turb % t_scale(c)          &
                    * ((  uu % n(c) * t % x(c)          &
                        + uv % n(c) * t % y(c)          &
                        + uw % n(c) * t % z(c))         &
                + afm_eta * (  ut % n(c) * u % x(c)     &
                             + vt % n(c) * u % y(c)     &
                             + wt % n(c) * u % z(c))    &
                + afm_psi * Flow % beta * Flow % grav_x * t2 % n(c))

      vt % n(c) = -c_theta * Turb % t_scale(c)          &
                    * ((  uv % n(c) * t % x(c)          &
                        + vv % n(c) * t % y(c)          &
                        + vw % n(c) * t % z(c))         &
                + afm_eta * (  ut % n(c) * v % x(c)     &
                             + vt % n(c) * v % y(c)     &
                             + wt % n(c) * v % z(c))    &
                + afm_psi * Flow % beta * Flow % grav_y * t2 % n(c))

      wt % n(c) = -c_theta * Turb % t_scale(c)          &
                    * ((  uw % n(c) * t % x(c)          &
                        + vw % n(c) * t % y(c)          &
                        + ww % n(c) * t % z(c))         &
                + afm_eta * (  ut % n(c) * w % x(c)     &
                             + vt % n(c) * w % y(c)     &
                             + wt % n(c) * w % z(c))    &
                + afm_psi * Flow % beta * Flow % grav_z * t2 % n(c))
    end do
  end if

  if(Turb % model .eq. K_EPS        .or.  &
     Turb % model .eq. K_EPS_ZETA_F .or.  &
     Turb % model .eq. HYBRID_LES_RANS) then

    do s = 1, Grid % n_faces
      c1 = Grid % faces_c(1,s)
      c2 = Grid % faces_c(2,s)

      if(c2 < 0) then

        if(Grid % Bnd_Cond_Type(c2) .eq. WALL .or. &
           Grid % Bnd_Cond_Type(c2) .eq. WALLFL) then

          nx = Grid % sx(s) / Grid % s(s)
          ny = Grid % sy(s) / Grid % s(s)
          nz = Grid % sz(s) / Grid % s(s)
          qx = t % q(c2) * nx
          qy = t % q(c2) * ny
          qz = t % q(c2) * nz

          ebf = Turb % Ebf_Momentum(c1)

          ut_log_law = - Turb % con_w(c1)                             &
                    / (Flow % density(c1) * Flow % capacity(c1))      &
                    * (t % n(c2) - t % n(c1)) / Grid % wall_dist(c1)  &
                    * nx
          vt_log_law = - Turb % con_w(c1)                             &
                    / (Flow % density(c1) * Flow % capacity(c1))      &
                    * (t % n(c2) - t % n(c1)) / Grid % wall_dist(c1)  &
                    * ny
          wt_log_law = - Turb % con_w(c1)                             &
                    / (Flow % density(c1) * Flow % capacity(c1))      &
                    * (t % n(c2) - t % n(c1)) / Grid % wall_dist(c1)  &
                    * nz

          ut % n(c1) = ut % n(c1) * exp(-1.0 * ebf)  &
                     + ut_log_law * exp(-1.0 / ebf)
          vt % n(c1) = vt % n(c1) * exp(-1.0 * ebf)  &
                     + vt_log_law * exp(-1.0 / ebf)
          wt % n(c1) = wt % n(c1) * exp(-1.0 * ebf)  &
                     + wt_log_law * exp(-1.0 / ebf)
        end if
      end if
    end do
  end if

  end subroutine
