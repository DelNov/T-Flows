!==============================================================================!
  subroutine Grad_Variable_With_Front(Vof, var, phif)
!------------------------------------------------------------------------------!
!   Calculates gradient of a variable from field Flow                          !
!                                                                              !
!   (Closely related (derived from) to Field_Mod_Grad_Variable)                !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Vof_Type)   :: Vof
  type(Var_Type)    :: var
  real, intent(in)  :: phif
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: Grid
  type(Field_Type), pointer :: Flow
!==============================================================================!

  ! Take alias
  Grid => Vof % pnt_grid
  Flow => Vof % pnt_flow

  ! Refresh buffers for variable
  call Grid % Exchange_Cells_Real(var % n)

  ! Compute individual gradients without refreshing buffers
  call Grad_Component_No_Refresh_With_Front(Vof, var % n, 1, var % x, phif)
  call Grad_Component_No_Refresh_With_Front(Vof, var % n, 2, var % y, phif)
  call Grad_Component_No_Refresh_With_Front(Vof, var % n, 3, var % z, phif)

  ! Refresh buffers for gradient components
  call Grid % Exchange_Cells_Real(var % x)
  call Grid % Exchange_Cells_Real(var % y)
  call Grid % Exchange_Cells_Real(var % z)

  end subroutine
