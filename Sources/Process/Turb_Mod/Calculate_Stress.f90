!==============================================================================!
  subroutine Calculate_Stress(Turb)
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
  class(Turb_Type), target :: Turb
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
  Flow => Turb % pnt_flow
  Grid => Flow % pnt_grid
  nc = Grid % n_cells
  nb = Grid % n_bnd_cells
  call Flow % Alias_Momentum    (u, v, w)
  call Turb % Alias_Stresses    (uu, vv, ww, uv, uw, vw)
  call Turb % Alias_K_Eps_Zeta_F(kin, eps, zeta, f22)

  call Flow % Grad_Variable(u)
  call Flow % Grad_Variable(v)
  call Flow % Grad_Variable(w)

  if( Turb % model .eq. K_EPS ) then
    do c = 1, Grid % n_cells

      uu % n(c) = - 2. * Turb % vis_t(c) / Flow % density(c)  &
                       * u % x(c) + TWO_THIRDS * kin % n(c)
      vv % n(c) = - 2. * Turb % vis_t(c) / Flow % density(c)  &
                       * v % y(c) + TWO_THIRDS * kin % n(c)
      ww % n(c) = - 2. * Turb % vis_t(c) / Flow % density(c)  &
                       * w % z(c) + TWO_THIRDS * kin % n(c)

      uv % n(c) = - Turb % vis_t(c) / Flow % density(c) * (u % y(c) + v % x(c))
      uw % n(c) = - Turb % vis_t(c) / Flow % density(c) * (u % z(c) + w % x(c))
      vw % n(c) = - Turb % vis_t(c) / Flow % density(c) * (v % z(c) + w % y(c))

    end do
  end if

  if(Turb % model .eq. K_EPS_ZETA_F .or.  &
     Turb % model .eq. HYBRID_LES_RANS) then

    call Flow % Grad(Grid % wall_dist, wd_x(-nb:nc),  &
                                       wd_y(-nb:nc),  &
                                       wd_z(-nb:nc))

    do c = 1, Grid % n_cells

      wd_m    = sqrt(wd_x(c)**2 + wd_y(c)**2 + wd_z(c)**2)
      wd_x(c) = abs(wd_x(c) / wd_m)
      wd_y(c) = abs(wd_y(c) / wd_m)
      wd_z(c) = abs(wd_z(c) / wd_m)

      ! Projections of v2 from k_eps_zeta_f model in three coordinate directions
      u2 = (zeta % n(c) * kin % n(c)) * wd_x(c)
      v2 = (zeta % n(c) * kin % n(c)) * wd_y(c)
      w2 = (zeta % n(c) * kin % n(c)) * wd_z(c)

      ! Take the part from Boussinesq's hypothesis (as if in k_eps) ...
      if(Turb % model .eq. K_EPS_ZETA_F) then 
        uu % n(c) = - 2. * Turb % vis_t(c) / Flow % density(c)  &
                         * u % x(c) + TWO_THIRDS * kin % n(c)
        vv % n(c) = - 2. * Turb % vis_t(c) / Flow % density(c)  &
                         * v % y(c) + TWO_THIRDS * kin % n(c)
        ww % n(c) = - 2. * Turb % vis_t(c) / Flow % density(c)  &
                         * w % z(c) + TWO_THIRDS * kin % n(c)

        uv % n(c) = - Turb % vis_t(c) / Flow % density(c) * (u % y(c) + v % x(c))
        uw % n(c) = - Turb % vis_t(c) / Flow % density(c) * (u % z(c) + w % x(c))
        vw % n(c) = - Turb % vis_t(c) / Flow % density(c) * (v % z(c) + w % y(c))
      else if(Turb % model .eq. HYBRID_LES_RANS) then
        uu % n(c) = - 2. * Turb % vis_t_eff(c) / Flow % density(c)  &
                         * u % x(c) + TWO_THIRDS * kin % n(c)
        vv % n(c) = - 2. * Turb % vis_t_eff(c) / Flow % density(c)  &
                         * v % y(c) + TWO_THIRDS * kin % n(c)
        ww % n(c) = - 2. * Turb % vis_t_eff(c) / Flow % density(c)  &
                         * w % z(c) + TWO_THIRDS * kin % n(c)

        uv % n(c) = - Turb % vis_t_eff(c) / Flow % density(c) * (u % y(c) + v % x(c))
        uw % n(c) = - Turb % vis_t_eff(c) / Flow % density(c) * (u % z(c) + w % x(c))
        vw % n(c) = - Turb % vis_t_eff(c) / Flow % density(c) * (v % z(c) + w % y(c))
      end if
 
      ! ... and balance them with explicitly resolved normal stresses
      uu % n(c) = uu % n(c) * (1.0 - wd_x(c)) + u2
      vv % n(c) = vv % n(c) * (1.0 - wd_y(c)) + v2
      ww % n(c) = ww % n(c) * (1.0 - wd_z(c)) + w2

    end do
  end if

  end subroutine
