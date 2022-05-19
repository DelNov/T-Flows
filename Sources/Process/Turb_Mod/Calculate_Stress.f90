!==============================================================================!
  subroutine Turb_Mod_Calculate_Stress(turb)
!------------------------------------------------------------------------------!
!   Calculates algebraic Reynolds stresses                                     !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Work_Mod, only: wd_x => r_cell_01,  &
                      wd_y => r_cell_02,  &
                      wd_z => r_cell_03
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Turb_Type),  target :: turb
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type), pointer :: Flow
  type(Grid_Type),  pointer :: Grid
  type(Var_Type),   pointer :: u, v, w
  type(Var_Type),   pointer :: uu, vv, ww, uv, uw, vw
  type(Var_Type),   pointer :: kin, eps, zeta, f22
  integer                   :: c, nc, nb
  real                      :: wd_m, u2, v2, w2
!==============================================================================!

  ! Take aliases
  Flow => turb % pnt_flow
  Grid => Flow % pnt_grid
  nc = Grid % n_cells
  nb = Grid % n_bnd_cells
  call Flow % Alias_Momentum(u, v, w)
  call Turb_Mod_Alias_Stresses    (turb, uu, vv, ww, uv, uw, vw)
  call Turb_Mod_Alias_K_Eps_Zeta_F(turb, kin, eps, zeta, f22)

  call Flow % Grad_Variable(u)
  call Flow % Grad_Variable(v)
  call Flow % Grad_Variable(w)

  if( turb % model .eq. K_EPS        .or.  &
      turb % model .eq. K_EPS_ZETA_F) then 

    do c = 1, Grid % n_cells

      uu % n(c) = - 2. * turb % vis_t(c) / Flow % density(c)  &
                       * u % x(c) + TWO_THIRDS * kin % n(c)
      vv % n(c) = - 2. * turb % vis_t(c) / Flow % density(c)  &
                       * v % y(c) + TWO_THIRDS * kin % n(c)
      ww % n(c) = - 2. * turb % vis_t(c) / Flow % density(c)  &
                       * w % z(c) + TWO_THIRDS * kin % n(c)

      uv % n(c) = - turb % vis_t(c) / Flow % density(c) * (u % y(c) + v % x(c))
      uw % n(c) = - turb % vis_t(c) / Flow % density(c) * (u % z(c) + w % x(c))
      vw % n(c) = - turb % vis_t(c) / Flow % density(c) * (v % z(c) + w % y(c))

    end do

  else if(turb % model .eq. HYBRID_LES_RANS) then

    do c = 1, Grid % n_cells

      uu % n(c) = - 2. * turb % vis_t_eff(c) / Flow % density(c)  &
                       * u % x(c) + TWO_THIRDS * kin % n(c)
      vv % n(c) = - 2. * turb % vis_t_eff(c) / Flow % density(c)  &
                       * v % y(c) + TWO_THIRDS * kin % n(c)
      ww % n(c) = - 2. * turb % vis_t_eff(c) / Flow % density(c)  &
                       * w % z(c) + TWO_THIRDS * kin % n(c)

      uv % n(c) = - turb % vis_t_eff(c) / Flow % density(c) * (u % y(c) + v % x(c))
      uw % n(c) = - turb % vis_t_eff(c) / Flow % density(c) * (u % z(c) + w % x(c))
      vw % n(c) = - turb % vis_t_eff(c) / Flow % density(c) * (v % z(c) + w % y(c))
    end do
  end if

  if( turb % model .eq. K_EPS_ZETA_F .or.  &
      turb % model .eq. HYBRID_LES_RANS) then
    do c = 1, Grid % n_cells

      uu % n(c) = zeta % n(c) * kin % n(c) 
      vv % n(c) = zeta % n(c) * kin % n(c)
      ww % n(c) = zeta % n(c) * kin % n(c)

    end do
  end if

  end subroutine
