!==============================================================================!
  subroutine Turb_Mod_Vis_T_Spalart_Allmaras(turb)
!------------------------------------------------------------------------------!
!   Computes the turbulent viscosity for RANS models.                          !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Turb_Type), target :: turb
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type), pointer :: Flow
  type(Grid_Type),  pointer :: grid
  type(Var_Type),   pointer :: u, v, w
  type(Var_Type),   pointer :: vis
  integer                   :: c
  real                      :: x_rat, f_v1
!==============================================================================!

  ! Take aliases
  Flow => turb % pnt_flow
  grid => Flow % pnt_grid
  vis  => turb % vis
  call Flow % Alias_Momentum(u, v, w)

  if(turb % model .eq. DES_SPALART) then
    do c = 1, grid % n_cells
      x_rat    = vis % n(c) / Flow % viscosity(c)
      f_v1     = x_rat**3/(x_rat**3 + c_v1**3)
      turb % vis_t(c) = Flow % density(c) * f_v1 * vis % n(c)
    end do
  end if

  if(turb % model .eq. SPALART_ALLMARAS) then
    do c = 1, grid % n_cells
      x_rat = vis % n(c) / Flow % viscosity(c)
      f_v1  = x_rat**3/(x_rat**3 + c_v1**3)
      turb % vis_t(c) = Flow % density(c) * f_v1 * vis % n(c)
    end do
  end if

  call Grid_Mod_Exchange_Cells_Real(grid, turb % vis_t)

  end subroutine
