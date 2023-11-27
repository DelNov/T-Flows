!==============================================================================!
  subroutine Destroy_Petsc(Pet, var_name)
!------------------------------------------------------------------------------!
!>  The subroutine Destroy_Petsc in the Petsc_Mod module is responsible for
!>  properly deallocating and destroying PETSc objects associated with a
!>  given variable in T-Flows.
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * Nullify pointers:                                                        !
!     - Sets the pointer to the grid (pnt_grid) in the Pet object to null.     !
!       This action disconnects the PETSc object from its associated grid.     !
!   * Resetting variables:                                                     !
!     - Resets the total number of unknowns (m_upper) and the number of        !
!       unknowns in the current processor (m_lower) to zero.                   !
!   * Destroying PETSc objects:                                                !
!     - Destroys the PETSc matrix (A)                                          !
!     - Destroys the PETSc vectors (x and b)                                   !
!     - Destroys the PETSc KSP solver context (ksp)                            !
!   * Logging and process indication:                                          !
!     - If running on the first process, it prints messages indicating the     !
!       start and completion of the PETSc destruction process for a specific   !
!       variable, given by var_name.                                           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Petsc_Type) :: Pet       !! parent object of the Petsc_Type
  character(VL)     :: var_name  !! variable name
!==============================================================================!

  nullify(Pet % pnt_grid)

  if(First_Proc()) then
    write(*,'(a,a4,a1)', advance='no')  &
      ' # Destroying PETSc for ', trim(var_name), ' '
  end if

  ! Nullify otal number of unknowns and unknowns in this processor only
  Pet % m_upper = 0
  Pet % m_lower = 0

  !---------------------------!
  !    Destroy PETSc matrix   !
  !---------------------------!
  call C_Petsc_Mat_Destroy(Pet % A)

  !---------------------------!
  !   Destroy PETSc vectors   !
  !---------------------------!
  call C_Petsc_Vec_Destroy(Pet % x)
  call C_Petsc_Vec_Destroy(Pet % b)

  !--------------------------!
  !   Destroy PETSc solver   !
  !--------------------------!
  call C_Petsc_Ksp_Destroy(Pet % ksp)

  if(First_Proc()) print *, 'done!'

  end subroutine

