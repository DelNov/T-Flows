!==============================================================================!
  subroutine Update_Physical_Properties(Vof)
!------------------------------------------------------------------------------!
!   Update physical properties based on volume fraction variable               !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Vof_Type), target :: Vof
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: Grid
  type(Field_Type), pointer :: Flow
  type(Var_Type),   pointer :: col
  integer                   :: c
!==============================================================================!

  ! Take aliases
  Flow => Vof  % pnt_flow
  Grid => Flow % pnt_grid
  col  => Vof % fun
  ! col  => Vof % smooth

  ! Density and viscosity in cells
  do c = -Grid % n_bnd_cells, Grid % n_cells
    Flow % density(c)      = col % n(c)         * Vof % phase_dens(1)  &
                           + (1.0 - col % n(c)) * Vof % phase_dens(2)
    Flow % viscosity(c)    = col % n(c)         * Vof % phase_visc(1)  &
                           + (1.0 - col % n(c)) * Vof % phase_visc(2)
    Flow % capacity(c)     = col % n(c)         * Vof % phase_capa(1)  &
                           + (1.0 - col % n(c)) * Vof % phase_capa(2)
    Flow % conductivity(c) = col % n(c)         * Vof % phase_cond(1)  &
                           + (1.0 - col % n(c)) * Vof % phase_cond(2)
  end do

  call Grid % Exchange_Cells_Real(Flow % density)
  call Grid % Exchange_Cells_Real(Flow % viscosity)
  call Grid % Exchange_Cells_Real(Flow % capacity)
  call Grid % Exchange_Cells_Real(Flow % conductivity)

  end subroutine
