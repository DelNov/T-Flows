!==============================================================================!
  subroutine Grad(Flow, Grid, phi, phi_x, phi_y, phi_z)
!------------------------------------------------------------------------------!
!>  Calculates gradient of generic variable phi by the least squares method
!>  and stores its three components in arrays phi_x, phi_y and phi_z.  This
!>  subroutine refreshes the buffers of the variable before calculating the
!>  gradients, and of the calculated gradient components.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Field_Type), intent(in)    :: Flow  !! parent flow object
  type(Grid_Type),   intent(in)    :: Grid  !! grid object
  real,              intent(inout) :: phi  (-Grid % n_bnd_cells:Grid % n_cells)
    !! field whose gradients are being calculated
  real,              intent(out)   :: phi_x(-Grid % n_bnd_cells:Grid % n_cells)
    !! x component of the calculated gradient
  real,              intent(out)   :: phi_y(-Grid % n_bnd_cells:Grid % n_cells)
    !! y component of the calculated gradient
  real,              intent(out)   :: phi_z(-Grid % n_bnd_cells:Grid % n_cells)
    !! z component of the calculated gradient
!==============================================================================!

  ! Refresh buffers for array
  call Grid % Exchange_Cells_Real(phi)

  ! Compute individual gradients without refreshing buffers
  call Flow % Grad_Component_No_Refresh(Grid, phi, 1, phi_x)  ! dp/dx
  call Flow % Grad_Component_No_Refresh(Grid, phi, 2, phi_y)  ! dp/dy
  call Flow % Grad_Component_No_Refresh(Grid, phi, 3, phi_z)  ! dp/dz

  ! Refresh buffers for gradient components
  call Grid % Exchange_Cells_Real(phi_x)
  call Grid % Exchange_Cells_Real(phi_y)
  call Grid % Exchange_Cells_Real(phi_z)

  end subroutine
