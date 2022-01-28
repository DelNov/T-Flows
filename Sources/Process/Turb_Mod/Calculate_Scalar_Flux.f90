!==============================================================================!
  subroutine Turb_Mod_Calculate_Scalar_Flux(turb, sc)
!------------------------------------------------------------------------------!
!   Computes turbulent heat fluxes                                             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Turb_Type),  target :: turb
  integer, intent(in)      :: sc
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type), pointer :: Flow
  type(Grid_Type),  pointer :: Grid
  type(Var_Type),   pointer :: uu, vv, ww, uv, uw, vw, kin, zeta, eps, f
  type(Var_Type),   pointer :: u, v, w
  type(Var_Type),   pointer :: phi
  integer                   :: c, k, c1, c2, s
  real                      :: uc_log_law, vc_log_law, wc_log_law
  real                      :: nx, ny, nz, qx, qy, qz, ebf
  real                      :: uc_new, vc_new, wc_new
!==============================================================================!

  ! Take aliases
  Flow => turb % pnt_flow
  Grid => Flow % pnt_grid
  phi  => Flow % scalar(sc)
  call Turb_Mod_Alias_Stresses   (turb, uu, vv, ww, uv, uw, vw)
  call Flow % Alias_Momentum     (u, v, w)
  call Turb_Mod_Alias_K_Eps_Zeta_F(turb, kin, eps, zeta, f)

  ! Check if these are already computed somewhere, ...
  ! ... maybe this call is not needed
  call Flow % Grad_Variable(phi)

  !-----------------------------------------!
  !   Compute the sources in the interior   !
  !-----------------------------------------!
  call Control_Mod_Turbulent_Schmidt_Number(sc_t)

  !-----------------------------------------!
  ! First guess is the flux defined by SGDH !
  !-----------------------------------------!
  do c = 1, Grid % n_cells

    turb % uc(c) = - turb % vis_t(c) / Flow % density(c) / sc_t * phi % x(c)
    turb % vc(c) = - turb % vis_t(c) / Flow % density(c) / sc_t * phi % y(c)
    turb % wc(c) = - turb % vis_t(c) / Flow % density(c) / sc_t * phi % z(c)

    if(turb % model .eq. HYBRID_LES_RANS) then
      turb % uc(c) = - turb % vis_t_eff(c) / Flow % density(c) &
                                        / sc_t * phi % x(c)
      turb % vc(c) = - turb % vis_t_eff(c) / Flow % density(c) &
                                        / sc_t * phi % y(c)
      turb % wc(c) = - turb % vis_t_eff(c) / Flow % density(c) &
                                        / sc_t * phi % z(c)
    end if
  end do


  if(turb % scalar_flux_model .eq. GGDH) then

    do c = 1, Grid % n_cells
      turb % uc(c) = -c_theta * turb % t_scale(c) * (uu % n(c) * phi % x(c)  +  &
                                                     uv % n(c) * phi % y(c)  +  &
                                                     uw % n(c) * phi % z(c))
      turb % vc(c) = -c_theta * turb % t_scale(c) * (uv % n(c) * phi % x(c)  +  &
                                                     vv % n(c) * phi % y(c)  +  &
                                                     vw % n(c) * phi % z(c))
      turb % wc(c) = -c_theta * turb % t_scale(c) * (uw % n(c) * phi % x(c)  +  &
                                                     vw % n(c) * phi % y(c)  +  &
                                                     ww % n(c) * phi % z(c))
    end do

  else if(turb % scalar_flux_model .eq. AFM) then
    call Flow % Grad_Variable(Flow % u)
    call Flow % Grad_Variable(Flow % v)
    call Flow % Grad_Variable(Flow % w)
    do k = 1, 3
      do c = 1, Grid % n_cells

        uc_new = -c_theta * turb % t_scale(c) * ((  uu % n(c) * phi % x(c)    &
                                                  + uv % n(c) * phi % y(c)    &
                                                  + uw % n(c) * phi % z(c))   &
                                        + 0.6*(  turb % uc(c) * u % x(c)      &
                                               + turb % vc(c) * u % y(c)      &
                                               + turb % wc(c) * u % z(c)))


        vc_new = -c_theta * turb % t_scale(c) * ((  uv % n(c) * phi % x(c)    &
                                                  + vv % n(c) * phi % y(c)    &
                                                  + vw % n(c) * phi % z(c))   &
                                        + 0.6*(  turb % uc(c) * v % x(c)      &
                                               + turb % vc(c) * v % y(c)      &
                                               + turb % wc(c) * v % z(c)))

        wc_new = -c_theta * turb % t_scale(c) * ((  uw % n(c) * phi % x(c)    &
                                                  + vw % n(c) * phi % y(c)    &
                                                  + ww % n(c) * phi % z(c))   &
                                        + 0.6*(  turb % uc(c) * w % x(c)      &
                                               + turb % vc(c) * w % y(c)      &
                                               + turb % wc(c) * w % z(c)))

        turb % uc(c) = turb % uc(c) * 0.7 + uc_new * 0.3
        turb % vc(c) = turb % vc(c) * 0.7 + vc_new * 0.3
        turb % wc(c) = turb % wc(c) * 0.7 + wc_new * 0.3

      end do
    end do
  end if

  if(turb % model .eq. K_EPS        .or.  &
     turb % model .eq. K_EPS_ZETA_F .or.  &
     turb % model .eq. HYBRID_LES_RANS) then

    do s = 1, Grid % n_faces
      c1 = Grid % faces_c(1,s)
      c2 = Grid % faces_c(2,s)

      if(c2 < 0) then

        if(Grid % Bnd_Cond_Type(c2) .eq. WALL .or. &
           Grid % Bnd_Cond_Type(c2) .eq. WALLFL) then

          nx = Grid % sx(s) / Grid % s(s)
          ny = Grid % sy(s) / Grid % s(s)
          nz = Grid % sz(s) / Grid % s(s)

          ebf = Turb_Mod_Ebf_Momentum(turb, c1)

          uc_log_law = - turb % diff_w(c1)  &
                    / Flow % density(c1)    &
                    * (phi % n(c2) - phi % n(c1)) / Grid % wall_dist(c1) * nx
          vc_log_law = - turb % diff_w(c1)  &
                    / Flow % density(c1)    &
                    * (phi % n(c2) - phi % n(c1)) / Grid % wall_dist(c1) * ny
          wc_log_law = - turb % diff_w(c1)  &
                    / Flow % density(c1)    &
                    * (phi % n(c2) - phi % n(c1)) / Grid % wall_dist(c1) * nz

          turb % uc(c1) = turb % uc(c1) * exp(-1.0 * ebf)  &
                           + uc_log_law * exp(-1.0 / ebf)
          turb % vc(c1) = turb % vc(c1) * exp(-1.0 * ebf)  &
                           + vc_log_law * exp(-1.0 / ebf)
          turb % wc(c1) = turb % wc(c1) * exp(-1.0 * ebf)  &
                           + wc_log_law * exp(-1.0 / ebf)
        end if
      end if
    end do
  end if

  end subroutine
