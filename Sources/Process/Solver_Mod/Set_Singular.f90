!==============================================================================!
  subroutine Set_Singular(Sol, phi)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Solver_Type) :: Sol
  type(Var_Type)     :: phi
!-----------------------------------[Locals]-----------------------------------!
  type(Matrix_Type), pointer :: A
  type(Grid_Type),   pointer :: Grid
  integer                    :: i
!==============================================================================!

  ! Take aliases
  A    => phi % pnt_matrix
  Grid => phi % pnt_grid

  ! Prepare solver for the singular matrix
  if(Sol % solvers == PETSC) then
    call C_Petsc_Mat_Set_Null_Space(phi % Pet % A)
  else
    do i = Cells_In_Domain()
      A % val(A % dia(i)) = A % val(A % dia(i)) * (1.0 + MICRO)
    end do
  end if

  end subroutine
