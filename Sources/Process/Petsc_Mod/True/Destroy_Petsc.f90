!==============================================================================!
  subroutine Destroy_Petsc(Pet, var_name)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Petsc_Type) :: Pet
  character(VL)     :: var_name
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

