!==============================================================================!
  subroutine Grad_Mod_Variable(var, boundary)
!------------------------------------------------------------------------------!
!   Calculates gradient of generic array.                                      !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Var_Type), target :: var
  logical                :: boundary
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: grid
  integer                  :: s, c, c1, c2, iter
!==============================================================================!

  ! Take aliases
  grid  => var % pnt_grid

  call Grad_Mod_Component(grid, var % n, 1, var % x, boundary)  ! dp/dx
  call Grad_Mod_Component(grid, var % n, 2, var % y, boundary)  ! dp/dy
  call Grad_Mod_Component(grid, var % n, 3, var % z, boundary)  ! dp/dz

  end subroutine
