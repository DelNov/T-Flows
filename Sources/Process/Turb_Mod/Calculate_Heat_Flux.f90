!==============================================================================!
  subroutine Turb_Mod_Calculate_Heat_Flux(turb)
!------------------------------------------------------------------------------!
!   Computes turbulent heat fluxes                                             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Turb_Type),  target :: turb
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type), pointer :: flow
  type(Grid_Type),  pointer :: grid
  type(Var_Type),   pointer :: uu, vv, ww, uv, uw, vw
  type(Var_Type),   pointer :: t, ut, vt, wt
  integer                   :: c
!==============================================================================!

  ! Take aliases
  flow => turb % pnt_flow
  grid => flow % pnt_grid
  t    => flow % t
  call Turb_Mod_Alias_Heat_Fluxes(turb, ut, vt, wt)
  call Turb_Mod_Alias_Stresses   (turb, uu, vv, ww, uv, uw, vw)

  ! Check if these are already computed somewhere, ...
  ! ... maybe this call is not needed
  call Field_Mod_Grad_Variable(flow, t)

  !-----------------------------------------!
  !   Compute the sources in the interior   !
  !-----------------------------------------!
  call Control_Mod_Turbulent_Prandtl_Number(pr_t)

  if(turb % heat_flux_model .eq. SGDH) then

    do c = 1, grid % n_cells
      pr_t = max(Turb_Mod_Prandtl_Number(turb, c), TINY)
      ut % n(c) = - turb % vis_t(c) / flow % density(c) / pr_t * t % x(c)
      vt % n(c) = - turb % vis_t(c) / flow % density(c) / pr_t * t % y(c)
      wt % n(c) = - turb % vis_t(c) / flow % density(c) / pr_t * t % z(c)

      if(turb % model .eq. HYBRID_LES_RANS) then
        ut % n(c) = - turb % vis_t_eff(c) / flow % density(c) &
                                          / pr_t * t % x(c)
        vt % n(c) = - turb % vis_t_eff(c) / flow % density(c) &
                                          / pr_t * t % y(c)
        wt % n(c) = - turb % vis_t_eff(c) / flow % density(c) &
                                          / pr_t * t % z(c)
      end if
    end do

  else if(turb % heat_flux_model .eq. GGDH) then

    do c = 1, grid % n_cells
      ut % n(c) = -c_theta * turb % t_scale(c) * (uu % n(c) * t % x(c)  +  &
                                                  uv % n(c) * t % y(c)  +  &
                                                  uw % n(c) * t % z(c))
      vt % n(c) = -c_theta * turb % t_scale(c) * (uv % n(c) * t % x(c)  +  &
                                                  vv % n(c) * t % y(c)  +  &
                                                  vw % n(c) * t % z(c))
      wt % n(c) = -c_theta * turb % t_scale(c) * (uw % n(c) * t % x(c)  +  &
                                                  vw % n(c) * t % y(c)  +  &
                                                  ww % n(c) * t % z(c))
    end do

  else if(turb % heat_flux_model .eq. AFM) then
  end if

  end subroutine
