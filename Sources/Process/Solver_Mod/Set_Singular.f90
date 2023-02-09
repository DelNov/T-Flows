!==============================================================================!
  subroutine Set_Singular(Sol, A)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Solver_Type) :: Sol
  type(Matrix_Type)  :: A
!-----------------------------------[Locals]-----------------------------------!
  integer :: i, ni
!==============================================================================!

  ! Prepare solver for the singular matrix
  if(Sol % solvers == PETSC) then
    call C_Petsc_Mat_Set_Null_Space(Sol % Pet % A)
  else
    ni = A % pnt_grid % n_cells - A % pnt_grid % Comm % n_buff_cells
    do i = 1, A % pnt_grid % n_cells
      A % val(A % dia(i)) = A % val(A % dia(i)) * (1.0 + MICRO)
    end do
  end if

  end subroutine
