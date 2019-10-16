!==============================================================================!
  subroutine Grad_Mod_Array(grid, phi, phi_x, phi_y, phi_z)
!------------------------------------------------------------------------------!
!   Calculates gradient of generic array.                                      !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  real            :: phi  (-grid % n_bnd_cells:grid % n_cells),  &
                     phi_x(-grid % n_bnd_cells:grid % n_cells),  &
                     phi_y(-grid % n_bnd_cells:grid % n_cells),  &
                     phi_z(-grid % n_bnd_cells:grid % n_cells)
!-----------------------------------[Locals]-----------------------------------!
  integer :: s, c, c1, c2, iter
!==============================================================================!

  call Grad_Mod_Component(grid, phi, 1, phi_x)  ! dphi/dx
  call Grad_Mod_Component(grid, phi, 2, phi_y)  ! dphi/dy
  call Grad_Mod_Component(grid, phi, 3, phi_z)  ! dphi/dz

  end subroutine
