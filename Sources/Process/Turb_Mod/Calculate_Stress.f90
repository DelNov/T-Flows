!==============================================================================!
  subroutine Turb_Mod_Calculate_Stress(turb)
!------------------------------------------------------------------------------!
!   Calculates algebraic Reynolds stresses                                     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Turb_Type),  target :: turb
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type), pointer :: flow
  type(Grid_Type),  pointer :: grid
  type(Var_Type),   pointer :: u, v, w
  type(Var_Type),   pointer :: uu, vv, ww, uv, uw, vw
  type(Var_Type),   pointer :: kin, eps
  integer                   :: c
!==============================================================================!

  ! Take aliases
  flow => turb % pnt_flow
  grid => flow % pnt_grid
  call Field_Mod_Alias_Momentum(flow, u, v, w)
  call Turb_Mod_Alias_Stresses (turb, uu, vv, ww, uv, uw, vw)
  call Turb_Mod_Alias_K_Eps    (turb, kin, eps)

  call Field_Mod_Grad_Variable(flow, u)
  call Field_Mod_Grad_Variable(flow, v)
  call Field_Mod_Grad_Variable(flow, w)

  do c = 1, grid % n_cells

    uu % n(c) = - 2. * turb % vis_t(c) * u % x(c) + TWO_THIRDS * kin % n(c)
    vv % n(c) = - 2. * turb % vis_t(c) * v % y(c) + TWO_THIRDS * kin % n(c)
    ww % n(c) = - 2. * turb % vis_t(c) * w % z(c) + TWO_THIRDS * kin % n(c)

    uv % n(c) = - turb % vis_t(c) * (u % y(c) + v % x(c))
    uw % n(c) = - turb % vis_t(c) * (u % z(c) + w % x(c))
    vw % n(c) = - turb % vis_t(c) * (v % z(c) + w % y(c))

  end do

  end subroutine
