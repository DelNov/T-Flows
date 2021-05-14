!==============================================================================!
  subroutine Multiphase_Mod_Vof_Grad_Variable_With_Front(mult, var, phif)
!------------------------------------------------------------------------------!
!   Calculates gradient of a variable from field flow                          !
!                                                                              !
!   (Closely related (derived from) to Field_Mod_Grad_Variable)                !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Multiphase_Type) :: mult
  type(Var_Type)        :: var
  real, intent(in)      :: phif
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: grid
  type(Field_Type), pointer :: flow
!==============================================================================!

  ! Take alias
  grid => mult % pnt_grid
  flow => mult % pnt_flow

  ! Refresh buffers for variable
  call Grid_Mod_Exchange_Cells_Real(grid, var % n)

  ! Compute individual gradients without refreshing buffers
  call Multiphase_Mod_Vof_Grad_Component_No_Refresh_With_Front(mult, var % n, 1, var % x, phif)
  call Multiphase_Mod_Vof_Grad_Component_No_Refresh_With_Front(mult, var % n, 2, var % y, phif)
  call Multiphase_Mod_Vof_Grad_Component_No_Refresh_With_Front(mult, var % n, 3, var % z, phif)

  ! Refresh buffers for gradient components
  call Grid_Mod_Exchange_Cells_Real(grid, var % x)
  call Grid_Mod_Exchange_Cells_Real(grid, var % y)
  call Grid_Mod_Exchange_Cells_Real(grid, var % z)

  end subroutine
