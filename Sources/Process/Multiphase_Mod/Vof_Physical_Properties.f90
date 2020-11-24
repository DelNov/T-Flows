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
  do c = 1, grid % n_cells
    flow % density(c)      = col % n(c)         * mult % phase_dens(1)  &
                           + (1.0 - col % n(c)) * mult % phase_dens(2)
    flow % viscosity(c)    = col % n(c)         * mult % phase_visc(1)  &
                           + (1.0 - col % n(c)) * mult % phase_visc(2)
    flow % capacity(c)     = col % n(c)         * mult % phase_capa(1)  &
                           + (1.0 - col % n(c)) * mult % phase_capa(2)
    flow % conductivity(c) = col % n(c)         * mult % phase_cond(1)  &
                           + (1.0 - col % n(c)) * mult % phase_cond(2)
  end do

  ! At boundaries (this shouldn't be needed with proper interpolation)
  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    if(c2 < 0) then
      if(Grid_Mod_Bnd_Cond_Type(grid, c2) .ne. INFLOW) then
        flow % density     (c2) = flow % density     (c1)
        flow % viscosity   (c2) = flow % viscosity   (c1)
        flow % capacity    (c2) = flow % capacity    (c1)
        flow % conductivity(c2) = flow % conductivity(c1)
      end if
    end if
  end do

  call Grid_Mod_Exchange_Cells_Real(grid, flow % density)
  call Grid_Mod_Exchange_Cells_Real(grid, flow % viscosity)
  call Grid_Mod_Exchange_Cells_Real(grid, flow % capacity)
  call Grid_Mod_Exchange_Cells_Real(grid, flow % conductivity)

  end subroutine
