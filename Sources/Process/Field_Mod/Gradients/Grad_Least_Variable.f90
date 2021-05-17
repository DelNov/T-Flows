!==============================================================================!
  subroutine Grad_Least_Variable(Flow, var)
!------------------------------------------------------------------------------!
!   Calculates gradient of a variable from field Flow with least squares       !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Field_Type) :: Flow
  type(Var_Type)    :: var
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: grid
!==============================================================================!

  ! Take alias
  grid => Flow % pnt_grid

  ! Refresh buffers for variable
  call Grid_Mod_Exchange_Cells_Real(grid, var % n)

  ! Compute individual gradients without refreshing buffers
  call Flow % Grad_Component_No_Refresh(var % n, 1, var % x)  ! dp/dx
  call Flow % Grad_Component_No_Refresh(var % n, 2, var % y)  ! dp/dy
  call Flow % Grad_Component_No_Refresh(var % n, 3, var % z)  ! dp/dz

  ! Refresh buffers for gradient components
  call Grid_Mod_Exchange_Cells_Real(grid, var % x)
  call Grid_Mod_Exchange_Cells_Real(grid, var % y)
  call Grid_Mod_Exchange_Cells_Real(grid, var % z)

  end subroutine
