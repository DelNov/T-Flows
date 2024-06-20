!==============================================================================!
  subroutine Add_Inertial_Term(Process, phi, Flow, Grid, coef_a, coef_b)
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
  class(Process_Type)      :: Process
  type(Var_Type),   target :: phi
  type(Field_Type), target :: Flow
  type(Grid_Type),  target :: Grid
  real                     :: coef_a(-Grid % n_bnd_cells:Grid % n_cells)
  real                     :: coef_b(-Grid % n_bnd_cells:Grid % n_cells)
!-----------------------------------[Locals]-----------------------------------!
  real, contiguous, pointer :: phi_o(:), b(:), vol(:)
  integer                   :: c
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Process)
!==============================================================================!

  call Profiler % Start('Add_Inertial_Term')

  ! Take some aliases
  b     => Flow % Nat % b
  phi_o => phi % o

  ! Unit for momentum: [kg/m^3 * 1 * m/s * m^3 / s = kg m / s^2 = N]
  ! Unit for energy:   [kg/m^3 * J/kg/K * K * m^3 / s = J / s = W]

  !$acc parallel loop independent                        &
  !$acc present(grid_region_f_cell, grid_region_l_cell,  &
  !$acc         grid_vol,                                &
  !$acc         b, phi_o, coef_a, coef_b)
  do c = Cells_In_Domain_Gpu()
    b(c) = b(c) + coef_a(c) * coef_b(c) * phi_o(c) * grid_vol(c) / Flow % dt
  end do
  !$acc end parallel

  call Profiler % Stop('Add_Inertial_Term')

  end subroutine
