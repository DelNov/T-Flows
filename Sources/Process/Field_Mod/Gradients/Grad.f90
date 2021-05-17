!==============================================================================!
  subroutine Grad(Flow, phi, phi_x, phi_y, phi_z)
!------------------------------------------------------------------------------!
!   Calculates gradient of an array defined in cells                           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Field_Type)  :: Flow
  real               :: phi  ( -Flow % pnt_grid % n_bnd_cells  &
                               :Flow % pnt_grid % n_cells)
  real               :: phi_x( -Flow % pnt_grid % n_bnd_cells  &
                               :Flow % pnt_grid % n_cells)
  real               :: phi_y( -Flow % pnt_grid % n_bnd_cells  &
                               :Flow % pnt_grid % n_cells)
  real               :: phi_z( -Flow % pnt_grid % n_bnd_cells  &
                               :Flow % pnt_grid % n_cells)
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: grid
!==============================================================================!

  ! Take alias
  grid => Flow % pnt_grid

  ! Refresh buffers for array
  call Grid_Mod_Exchange_Cells_Real(grid, phi)

  ! Compute individual gradients without refreshing buffers
  call Flow % Grad_Component_No_Refresh(phi, 1, phi_x)  ! dp/dx
  call Flow % Grad_Component_No_Refresh(phi, 2, phi_y)  ! dp/dy
  call Flow % Grad_Component_No_Refresh(phi, 3, phi_z)  ! dp/dz

  ! Refresh buffers for gradient components
  call Grid_Mod_Exchange_Cells_Real(grid, phi_x)
  call Grid_Mod_Exchange_Cells_Real(grid, phi_y)
  call Grid_Mod_Exchange_Cells_Real(grid, phi_z)

  end subroutine
