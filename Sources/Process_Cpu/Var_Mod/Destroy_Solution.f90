!==============================================================================!
  subroutine Var_Mod_Destroy_Solution(phi)
!------------------------------------------------------------------------------!
!>  This subroutine is responsible for deallocating and cleaning up a complex
!>  variable (created by Var_Mod_Create_Solution) whose value is obtained by
!>  solving partial differential equations.
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * Comprehensive cleanup: The subroutine methodically deallocates the       !
!     current, old, and older-than-old values, boundary values, fluxes,        !
!     gradient components, and boundary condition types.                       !
!   * PETSc integration: If PETSc is enabled, the subroutine also takes care   !
!     of destroying the PETSc instance associated with the variable. This is   !
!     crucial for maintaining the integrity of the PETSc environment and       !
!     freeing up resources used by PETSc.                                      !
!   * Finalization steps: After deallocating all data and destroying PETSc     !
!     instances (if applicable), the subroutine resets the variable's name     !
!     and flux name. This step signifies the complete removal of the variable  !
!     from the simulation.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Var_Type) :: phi  !! variable object being destroyed
!==============================================================================!

  ! Nullify pointers to grid and matrix
  nullify(phi % pnt_grid)
  nullify(phi % pnt_matrix)

  ! Values new (n), old (o), and older than old (oo)
  deallocate(phi % n)
  deallocate(phi % o)
  deallocate(phi % oo)

  ! Variable's boundary value
  deallocate(phi % b)

  ! Variable's boundary flux
  deallocate(phi % q)

  ! Boundary cell type
  deallocate(phi % bnd_cond_type)

  ! Gradients
  deallocate(phi % x)
  deallocate(phi % y)
  deallocate(phi % z)

  ! Destroy PETSc too
# if T_FLOWS_PETSC == 1
    call phi % Pet % Destroy_Petsc(phi % name)
# endif

  ! Delete variable name (do this last because name is used in Destroy_Petsc)
  phi % name      = ''
  phi % flux_name = ''

  end subroutine
