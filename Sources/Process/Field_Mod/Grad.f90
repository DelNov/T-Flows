!==============================================================================!
  subroutine Field_Mod_Grad(flow, phi, impose_symmetry)
!------------------------------------------------------------------------------!
!   Calculates gradient of an array defined in cells                           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type)  :: flow
  real              :: phi  (- flow % pnt_grid % n_bnd_cells  &
                             : flow % pnt_grid % n_cells)
  real              :: phi_x(- flow % pnt_grid % n_bnd_cells  &
                             : flow % pnt_grid % n_cells)
  real              :: phi_y(- flow % pnt_grid % n_bnd_cells  &
                            : flow % pnt_grid % n_cells)
  real              :: phi_z(- flow % pnt_grid % n_bnd_cells  &
                             : flow % pnt_grid % n_cells)
  logical, optional :: impose_symmetry
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: grid
!==============================================================================!

  ! Take alias
  grid => flow % pnt_grid

  ! Refresh buffers for array
  call Grid_Mod_Exchange_Cells_Real(grid, phi)

  ! Compute individual gradients without refreshing buffers
  if(present(impose_symmetry)) then
    call Field_Mod_Grad_Component_No_Refresh(flow, phi, 1, phi_x,  &
                                             impose_symmetry)  ! dp/dx
    call Field_Mod_Grad_Component_No_Refresh(flow, phi, 2, phi_y,  &
                                             impose_symmetry)  ! dp/dy
    call Field_Mod_Grad_Component_No_Refresh(flow, phi, 3, phi_z,  &
                                             impose_symmetry)  ! dp/dz
  else
    call Field_Mod_Grad_Component_No_Refresh(flow, phi, 1, phi_x)  ! dp/dx
    call Field_Mod_Grad_Component_No_Refresh(flow, phi, 2, phi_y)  ! dp/dy
    call Field_Mod_Grad_Component_No_Refresh(flow, phi, 3, phi_z)  ! dp/dz
  end if

  ! Refresh buffers for gradient components
  call Grid_Mod_Exchange_Cells_Real(grid, phi_x)
  call Grid_Mod_Exchange_Cells_Real(grid, phi_y)
  call Grid_Mod_Exchange_Cells_Real(grid, phi_z)

  end subroutine
