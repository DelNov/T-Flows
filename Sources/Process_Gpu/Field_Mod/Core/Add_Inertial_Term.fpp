!==============================================================================!
  subroutine Add_Inertial_Term(Flow, Grid, phi, coef)
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
  class(Field_Type), target :: Flow
  type(Grid_Type)           :: Grid
  type(Var_Type),    target :: phi
  real                      :: coef(-Grid % n_bnd_cells:Grid % n_cells)
!-----------------------------------[Locals]-----------------------------------!
  real, contiguous, pointer :: phi_o(:), b(:)
  integer                   :: c
!==============================================================================!

  call Profiler % Start('Add_Inertial_Term')

  ! Take some aliases
  b     => Flow % Nat % b
  phi_o => phi % o

  ! Unit for momentum: [kg/m^3 * 1 * m/s * m^3 / s = kg m / s^2 = N]
  ! Unit for energy:   [kg/m^3 * J/kg/K * K * m^3 / s = J / s = W]

  !$tf-acc loop begin
  do c = Cells_In_Domain()
    b(c) = b(c) + coef(c) * phi_o(c) * Grid % vol(c) / Flow % dt
  end do
  !$tf-acc loop end

  call Profiler % Stop('Add_Inertial_Term')

  end subroutine
