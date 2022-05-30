!==============================================================================!
  subroutine Calculate_Scalar_Flux(Turb, sc)
!------------------------------------------------------------------------------!
!   Computes turbulent heat fluxes                                             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Turb_Type), target :: Turb
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
  Flow => Turb % pnt_flow
  Grid => Flow % pnt_grid
  phi  => Flow % scalar(sc)
  call Turb % Alias_Stresses    (uu, vv, ww, uv, uw, vw)
  call Flow % Alias_Momentum    (u, v, w)
  call Turb % Alias_K_Eps_Zeta_F(kin, eps, zeta, f)

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
  if(Turb % scalar_flux_model .eq. SGDH) then
    do c = 1, Grid % n_cells

      Turb % uc(c) = - Turb % vis_t(c) / Flow % density(c) / sc_t * phi % x(c)
      Turb % vc(c) = - Turb % vis_t(c) / Flow % density(c) / sc_t * phi % y(c)
      Turb % wc(c) = - Turb % vis_t(c) / Flow % density(c) / sc_t * phi % z(c)

      if(Turb % model .eq. HYBRID_LES_RANS) then
        Turb % uc(c) = - Turb % vis_t_eff(c) / Flow % density(c) &
                                             / sc_t * phi % x(c)
        Turb % vc(c) = - Turb % vis_t_eff(c) / Flow % density(c) &
                                             / sc_t * phi % y(c)
        Turb % wc(c) = - Turb % vis_t_eff(c) / Flow % density(c) &
                                             / sc_t * phi % z(c)
      end if
    end do

  else if(Turb % scalar_flux_model .eq. GGDH) then

    do c = 1, Grid % n_cells
      Turb % uc(c) = -c_theta * Turb % t_scale(c) * (uu % n(c) * phi % x(c)  +  &
                                                     uv % n(c) * phi % y(c)  +  &
                                                     uw % n(c) * phi % z(c))
      Turb % vc(c) = -c_theta * Turb % t_scale(c) * (uv % n(c) * phi % x(c)  +  &
                                                     vv % n(c) * phi % y(c)  +  &
                                                     vw % n(c) * phi % z(c))
      Turb % wc(c) = -c_theta * Turb % t_scale(c) * (uw % n(c) * phi % x(c)  +  &
                                                     vw % n(c) * phi % y(c)  +  &
                                                     ww % n(c) * phi % z(c))
    end do

  else if(Turb % scalar_flux_model .eq. AFM) then
    call Flow % Grad_Variable(Flow % u)
    call Flow % Grad_Variable(Flow % v)
    call Flow % Grad_Variable(Flow % w)
    do k = 1, 3
      do c = 1, Grid % n_cells

        Turb % uc(c) = -c_theta*Turb % t_scale(c) * (( uu % n(c) * phi % x(c)    &
                                                     + uv % n(c) * phi % y(c)    &
                                                     + uw % n(c) * phi % z(c))   &
                                     + afm_eta * (  Turb % uc(c) * u % x(c)      &
                                                  + Turb % vc(c) * u % y(c)      &
                                                  + Turb % wc(c) * u % z(c)))


        Turb % vc(c) = -c_theta*Turb % t_scale(c) * (( uv % n(c) * phi % x(c)    &
                                                     + vv % n(c) * phi % y(c)    &
                                                     + vw % n(c) * phi % z(c))   &
                                     + afm_eta * (  Turb % uc(c) * v % x(c)      &
                                                  + Turb % vc(c) * v % y(c)      &
                                                  + Turb % wc(c) * v % z(c)))

        Turb % wc(c) = -c_theta*Turb % t_scale(c) * (( uw % n(c) * phi % x(c)    &
                                                     + vw % n(c) * phi % y(c)    &
                                                     + ww % n(c) * phi % z(c))   &
                                     + afm_eta * (  Turb % uc(c) * w % x(c)      &
                                                  + Turb % vc(c) * w % y(c)      &
                                                  + Turb % wc(c) * w % z(c)))

      end do
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

          ebf = Turb % Ebf_Momentum(c1)

          uc_log_law = - Turb % diff_w(c1)  &
                    / Flow % density(c1)    &
                    * (phi % n(c2) - phi % n(c1)) / Grid % wall_dist(c1) * nx
          vc_log_law = - Turb % diff_w(c1)  &
                    / Flow % density(c1)    &
                    * (phi % n(c2) - phi % n(c1)) / Grid % wall_dist(c1) * ny
          wc_log_law = - Turb % diff_w(c1)  &
                    / Flow % density(c1)    &
                    * (phi % n(c2) - phi % n(c1)) / Grid % wall_dist(c1) * nz

          Turb % uc(c1) = Turb % uc(c1) * exp(-1.0 * ebf)  &
                           + uc_log_law * exp(-1.0 / ebf)
          Turb % vc(c1) = Turb % vc(c1) * exp(-1.0 * ebf)  &
                           + vc_log_law * exp(-1.0 / ebf)
          Turb % wc(c1) = Turb % wc(c1) * exp(-1.0 * ebf)  &
                           + wc_log_law * exp(-1.0 / ebf)
        end if
      end if
    end do
  end if

  end subroutine
