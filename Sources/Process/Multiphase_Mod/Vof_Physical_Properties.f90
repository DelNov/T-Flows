!==============================================================================!
  subroutine Multiphase_Mod_Vof_Physical_Properties(mult)
!------------------------------------------------------------------------------!
!   Update physical properties based on volume fraction variable               !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Multiphase_Type), target :: mult
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: grid
  type(Field_Type), pointer :: flow
  type(Var_Type),   pointer :: col
  integer                   :: c, s, c1, c2
!==============================================================================!

  ! Take aliases
  flow  => mult % pnt_flow
  grid  => flow % pnt_grid
  col   => mult % vof
  ! col   => mult % smooth

  ! Density and viscosity in cells
  do c = -grid % n_bnd_cells, grid % n_cells
    flow % density(c)      = col % n(c)         * mult % phase_dens(1)  &
                           + (1.0 - col % n(c)) * mult % phase_dens(2)
    flow % viscosity(c)    = col % n(c)         * mult % phase_visc(1)  &
                           + (1.0 - col % n(c)) * mult % phase_visc(2)
    flow % capacity(c)     = col % n(c)         * mult % phase_capa(1)  &
                           + (1.0 - col % n(c)) * mult % phase_capa(2)
    flow % conductivity(c) = col % n(c)         * mult % phase_cond(1)  &
                           + (1.0 - col % n(c)) * mult % phase_cond(2)
  end do

  call Grid_Mod_Exchange_Cells_Real(grid, flow % density)
  call Grid_Mod_Exchange_Cells_Real(grid, flow % viscosity)
  call Grid_Mod_Exchange_Cells_Real(grid, flow % capacity)
  call Grid_Mod_Exchange_Cells_Real(grid, flow % conductivity)

  end subroutine
