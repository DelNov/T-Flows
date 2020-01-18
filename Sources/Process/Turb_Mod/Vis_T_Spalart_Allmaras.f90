!==============================================================================!
  subroutine Turb_Mod_Vis_T_Spalart_Allmaras(turb)
!------------------------------------------------------------------------------!
!   Computes the turbulent viscosity for RANS models.                          !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Turb_Type), target :: turb
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type), pointer :: flow
  type(Grid_Type),  pointer :: grid
  type(Var_Type),   pointer :: u, v, w
  type(Var_Type),   pointer :: vis
  integer                   :: c
  real                      :: x_rat, f_v1
!==============================================================================!

  ! Take aliases
  flow => turb % pnt_flow
  grid => flow % pnt_grid
  vis  => turb % vis
  call Field_Mod_Alias_Momentum(flow, u, v, w)

  if(turbulence_model .eq. DES_SPALART) then
    do c = 1, grid % n_cells
      x_rat    = vis % n(c) / flow % viscosity(c)
      f_v1     = x_rat**3/(x_rat**3 + c_v1**3)
      turb % vis_t(c) = flow % density(c) * f_v1 * vis % n(c)
    end do
  end if

  if(turbulence_model .eq. SPALART_ALLMARAS) then
    do c = 1, grid % n_cells
      x_rat = vis % n(c) / flow % viscosity(c)
      f_v1  = x_rat**3/(x_rat**3 + c_v1**3)
      turb % vis_t(c) = flow % density(c) * f_v1 * vis % n(c)
    end do
  end if

  call Grid_Mod_Exchange_Real(grid, turb % vis_t)

  end subroutine
