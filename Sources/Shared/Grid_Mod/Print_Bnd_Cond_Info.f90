!==============================================================================!
  subroutine Grid_Mod_Print_Bnd_Cond_List(grid)
!------------------------------------------------------------------------------!
!   Prints a list of boundary conditions in a grid.                            !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer :: j
!==============================================================================!

  print "(a)", " #======================================================"
  print "(a)", " # Grid currently has the following boundary conditions:"
  print "(a)", " #------------------------------------------------------"
  do j = 1, grid % n_bnd_cond
    print '(a3,i2,a2,a)', ' # ', j, '. ', grid % bnd_cond % name(j)
  end do
  print "(a)", " #------------------------------------------------------"

  end subroutine
