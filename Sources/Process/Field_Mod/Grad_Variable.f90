!==============================================================================!
  subroutine Field_Mod_Grad_Variable(flow, var)
!------------------------------------------------------------------------------!
!   Calculates gradient of a variable from field flow                          !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type) :: flow
  type(Var_Type)   :: var
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: grid
!==============================================================================!

  ! Take alias
  grid => flow % pnt_grid

  ! Refresh buffers for variable
  call Grid_Mod_Exchange_Cells_Real(grid, var % n)

  ! Compute individual gradients without refreshing buffers
  call Field_Mod_Grad_Component_No_Refresh(flow, var % n, 1, var % x)  ! dp/dx
  call Field_Mod_Grad_Component_No_Refresh(flow, var % n, 2, var % y)  ! dp/dy
  call Field_Mod_Grad_Component_No_Refresh(flow, var % n, 3, var % z)  ! dp/dz

  ! Refresh buffers for gradient components
  call Grid_Mod_Exchange_Cells_Real(grid, var % x)
  call Grid_Mod_Exchange_Cells_Real(grid, var % y)
  call Grid_Mod_Exchange_Cells_Real(grid, var % z)

  end subroutine
