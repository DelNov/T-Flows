!==============================================================================!
  subroutine Multiphase_Mod_Update_Physical_Properties(mult)
!------------------------------------------------------------------------------!
!    Update physical properties based on volume fraction variable              !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Multiphase_Type), target :: mult
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: grid
  type(Field_Type), pointer :: flow
  type(Var_Type),   pointer :: vof
  real,             pointer :: vof_f(:)
  integer                   :: c, s
!==============================================================================!

  ! Take aliases
  flow  => mult % pnt_flow
  grid  => flow % pnt_grid
  vof   => mult % vof
  vof_f => mult % vof_f

  ! Density and viscosity in cells
  do c = 1, grid % n_cells
    flow % density(c)   = vof % n(c)         * mult % phase_dens(1)      &
                        + (1.0 - vof % n(c)) * mult % phase_dens(2)
    flow % viscosity(c) = vof % n(c)         * mult % phase_visc(1)      &
                        + (1.0 - vof % n(c)) * mult % phase_visc(2)
  end do
  call Grid_Mod_Exchange_Real(grid, flow % density)
  call Grid_Mod_Exchange_Real(grid, flow % viscosity)

  ! Density at faces
  do s = 1, grid % n_faces
    flow % density_f(s) =        vof_f(s)  * mult % phase_dens(1)   &
                        + (1.0 - vof_f(s)) * mult % phase_dens(2)
  end do

  end subroutine
