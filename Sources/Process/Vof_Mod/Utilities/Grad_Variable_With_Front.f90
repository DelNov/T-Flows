!==============================================================================!
  subroutine Grad_Variable_With_Front(Vof, var, phif)
!------------------------------------------------------------------------------!
!   Calculates gradient of a variable from field flow                          !
!                                                                              !
!   (Closely related (derived from) to Field_Mod_Grad_Variable)                !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Vof_Type)   :: Vof
  type(Var_Type)    :: var
  real, intent(in)  :: phif
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: grid
  type(Field_Type), pointer :: flow
!==============================================================================!

  ! Take alias
  grid => Vof % pnt_grid
  flow => Vof % pnt_flow

  ! Refresh buffers for variable
  call Grid_Mod_Exchange_Cells_Real(grid, var % n)

  ! Compute individual gradients without refreshing buffers
  call Grad_Component_No_Refresh_With_Front(Vof, var % n, 1, var % x, phif)
  call Grad_Component_No_Refresh_With_Front(Vof, var % n, 2, var % y, phif)
  call Grad_Component_No_Refresh_With_Front(Vof, var % n, 3, var % z, phif)

  ! Refresh buffers for gradient components
  call Grid_Mod_Exchange_Cells_Real(grid, var % x)
  call Grid_Mod_Exchange_Cells_Real(grid, var % y)
  call Grid_Mod_Exchange_Cells_Real(grid, var % z)

  end subroutine
