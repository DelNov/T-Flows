!==============================================================================!
  subroutine Calculate_Stress(Turb)
!------------------------------------------------------------------------------!
!   Calculates algebraic Reynolds stresses                                     !
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

  if( Turb % model .eq. K_EPS        .or.  &
      Turb % model .eq. K_EPS_ZETA_F) then

    do c = Cells_In_Domain_And_Buffers()

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

  else if(Turb % model .eq. HYBRID_LES_RANS) then

    do c = Cells_In_Domain_And_Buffers()

      uu % n(c) = - 2. * Turb % vis_t_eff(c) / Flow % density(c)  &
                       * u % x(c) + TWO_THIRDS * kin % n(c)
      vv % n(c) = - 2. * Turb % vis_t_eff(c) / Flow % density(c)  &
                       * v % y(c) + TWO_THIRDS * kin % n(c)
      ww % n(c) = - 2. * Turb % vis_t_eff(c) / Flow % density(c)  &
                       * w % z(c) + TWO_THIRDS * kin % n(c)

      uv % n(c) = -Turb % vis_t_eff(c) / Flow % density(c) * (u % y(c)+v % x(c))
      uw % n(c) = -Turb % vis_t_eff(c) / Flow % density(c) * (u % z(c)+w % x(c))
      vw % n(c) = -Turb % vis_t_eff(c) / Flow % density(c) * (v % z(c)+w % y(c))
    end do
  end if

  if( Turb % model .eq. K_EPS_ZETA_F .or.  &
      Turb % model .eq. HYBRID_LES_RANS) then
    do c = Cells_In_Domain_And_Buffers()

      uu % n(c) = zeta % n(c) * kin % n(c)
      vv % n(c) = zeta % n(c) * kin % n(c)
      ww % n(c) = zeta % n(c) * kin % n(c)

    end do
  end if

  end subroutine
