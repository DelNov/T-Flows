!==============================================================================!
  subroutine Grad_Mod_Variable(var)
!------------------------------------------------------------------------------!
!   Calculates gradient of generic array.                                      !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Var_Type), target :: var
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: grid
  integer                  :: s, c, c1, c2, iter
!==============================================================================!

  ! Take aliases
  grid  => var % pnt_grid

  call Grad_Mod_Component(grid, var % n, 1, var % x)  ! dp/dx
  call Grad_Mod_Component(grid, var % n, 2, var % y)  ! dp/dy
  call Grad_Mod_Component(grid, var % n, 3, var % z)  ! dp/dz

  end subroutine
