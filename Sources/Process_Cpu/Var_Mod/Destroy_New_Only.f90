!==============================================================================!
  subroutine Var_Mod_Destroy_New_Only(phi)
!------------------------------------------------------------------------------!
!>  This subroutine is designed to deallocate and clean up a simplified
!>  variable (created by Var_Mod_Create_New_Only) at the end of its usage.
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * The subroutine deallocates the current value, gradients, and boundary    !
!     values of the variable. It also nullifies the pointer to the grid,       !
!     ensuring that the variable is completely disassociated from its          !
!     computational context.                                                   !
!   * Resetting Variable Attributes: The variable's name is reset, indicating  !
!     the removal of the variable from the simulation and its disconnection    !
!     from boundary and initial conditions defined in the control file.        !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Var_Type) :: phi  !! variable object being destroyed
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
