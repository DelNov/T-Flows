!==============================================================================!
  subroutine Set_Singular(Sol, phi)
!------------------------------------------------------------------------------!
!>  The Set_Singular subroutine in the Solver_Mod module plays an important
!>  role in preparing the solver for handling singular matrices. This
!>  functionality is especially relevant when dealing with matrices that might
!>  be deficient in rank or lack full invertibility, which is typically the
!>  case for the pressure correction solution.  Just like her sister "Run",
!>  the subroutine shows how encapsulation of decision making on wheather to
!>  use native or PETSc solvers works in practice in T-Flows.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Solver_Type) :: Sol  !! parent class
  type(Var_Type)     :: phi  !! unknown variable
!-----------------------------------[Locals]-----------------------------------!
  type(Matrix_Type), pointer :: A
  integer                    :: i, ni
!==============================================================================!

  ! Take aliases
  A => phi % pnt_matrix

  ! Prepare solver for the singular matrix
  if(Sol % solvers == PETSC) then
    call C_Petsc_Mat_Set_Null_Space(phi % Pet % A)
  else
    ni = A % pnt_grid % n_cells - A % pnt_grid % Comm % n_buff_cells
    do i = 1, A % pnt_grid % n_cells
      A % val(A % dia(i)) = A % val(A % dia(i)) * (1.0 + MICRO)
    end do
  end if

  end subroutine
