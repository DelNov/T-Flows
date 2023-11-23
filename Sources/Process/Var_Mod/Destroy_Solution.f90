!==============================================================================!
  subroutine Var_Mod_Destroy_Solution(phi)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Var_Type) :: phi
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
