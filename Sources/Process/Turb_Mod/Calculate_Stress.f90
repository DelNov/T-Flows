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
  type(Grid_Type),  pointer :: grid
  type(Var_Type),   pointer :: u, v, w
  type(Var_Type),   pointer :: uu, vv, ww, uv, uw, vw
  type(Var_Type),   pointer :: kin, eps, zeta, f22
  integer                   :: c, nc, nb
  real                      :: wd_m, u2, v2, w2
!==============================================================================!

  ! Take aliases
  Flow => turb % pnt_flow
  grid => Flow % pnt_grid
  nc = grid % n_cells
  nb = grid % n_bnd_cells
  call Flow % Alias_Momentum(u, v, w)
  call Turb_Mod_Alias_Stresses    (turb, uu, vv, ww, uv, uw, vw)
  call Turb_Mod_Alias_K_Eps_Zeta_F(turb, kin, eps, zeta, f22)

  call Flow % Grad_Variable(u)
  call Flow % Grad_Variable(v)
  call Flow % Grad_Variable(w)

  if( turb % model .eq. K_EPS ) then
    do c = 1, grid % n_cells

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
  end if

  if( turb % model .eq. K_EPS_ZETA_F ) then

    call Flow % Grad(grid % wall_dist, wd_x(-nb:nc),  &
                                       wd_y(-nb:nc),  &
                                       wd_z(-nb:nc))

    do c = 1, grid % n_cells

      wd_m = sqrt(wd_x(c)**2 + wd_y(c)**2 + wd_z(c)**2)
      wd_x(c) = abs(wd_x(c) / wd_m)
      wd_y(c) = abs(wd_y(c) / wd_m)
      wd_z(c) = abs(wd_z(c) / wd_m)

      ! Projections of v2 from k_eps_zeta_f model in three coordinate directions
      u2 = (zeta % n(c) * kin % n(c)) * wd_x(c)
      v2 = (zeta % n(c) * kin % n(c)) * wd_y(c)
      w2 = (zeta % n(c) * kin % n(c)) * wd_z(c)

      ! Take the part from Boussinesq's hypothesis (as if in k_eps) ...
      uu % n(c) = - 2. * turb % vis_t(c) / Flow % density(c)  &
                       * u % x(c) + TWO_THIRDS * kin % n(c)
      vv % n(c) = - 2. * turb % vis_t(c) / Flow % density(c)  &
                       * v % y(c) + TWO_THIRDS * kin % n(c)
      ww % n(c) = - 2. * turb % vis_t(c) / Flow % density(c)  &
                       * w % z(c) + TWO_THIRDS * kin % n(c)

      ! ... and balance them with explicitly resolved normal stresses
      uu % n(c) = uu % n(c) * (1.0 - wd_x(c)) + u2
      vv % n(c) = vv % n(c) * (1.0 - wd_y(c)) + v2
      ww % n(c) = ww % n(c) * (1.0 - wd_z(c)) + w2

      uv % n(c) = - turb % vis_t(c) / Flow % density(c) * (u % y(c) + v % x(c))
      uw % n(c) = - turb % vis_t(c) / Flow % density(c) * (u % z(c) + w % x(c))
      vw % n(c) = - turb % vis_t(c) / Flow % density(c) * (v % z(c) + w % y(c))

    end do
  end if

  end subroutine
