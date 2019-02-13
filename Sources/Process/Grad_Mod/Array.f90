!==============================================================================!
  subroutine Grad_Mod_Array(grid, phi, phi_x, phi_y, phi_z, boundary)
!------------------------------------------------------------------------------!
!   Calculates gradient of generic array.                                      !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Comm_Mod
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  real            :: phi  (-grid % n_bnd_cells:grid % n_cells),  &
                     phi_x(-grid % n_bnd_cells:grid % n_cells),  &
                     phi_y(-grid % n_bnd_cells:grid % n_cells),  &
                     phi_z(-grid % n_bnd_cells:grid % n_cells)
  logical         :: boundary
!-----------------------------------[Locals]-----------------------------------!
  integer :: s, c, c1, c2, iter
!==============================================================================!

  call Grad_Mod_Component(grid, phi, 1, phi_x, boundary)  ! dp/dx
  call Grad_Mod_Component(grid, phi, 2, phi_y, boundary)  ! dp/dy
  call Grad_Mod_Component(grid, phi, 3, phi_z, boundary)  ! dp/dz

  end subroutine
