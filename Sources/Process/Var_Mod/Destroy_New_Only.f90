!==============================================================================!
  subroutine Var_Mod_Destroy_New_Only(phi)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Var_Type) :: phi
!==============================================================================!

  ! Nullify the pointer to the grid
  nullify(phi % pnt_grid)

  ! Erase the variable name
  phi % name = ''

  ! Deallocate new value
  deallocate(phi % n)

  ! Deallocate gradients
  deallocate(phi % x)
  deallocate(phi % y)
  deallocate(phi % z)

  ! Deallocate variable's boundary value
  deallocate(phi % b)

  end subroutine
