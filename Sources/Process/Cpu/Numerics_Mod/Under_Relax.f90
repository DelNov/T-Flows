!==============================================================================!
  subroutine Numerics_Mod_Under_Relax(phi, A, b)
!------------------------------------------------------------------------------!
!>  This subroutine is a part of the Numerics_Mod module in T-Flows, designed
!>  to implement under-relaxation on a system of equations prior to solving
!>  with a linear solver.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Var_Type)    :: phi   !! variable to be solved
  type(Matrix_Type) :: A     !! linear system matrix
  real              :: b(:)  !! right-hand side vector
!-----------------------------------[Locals]-----------------------------------!
  integer :: c
!==============================================================================!

  !------------------------------------------!
  !   Browse through cells and under-relax   !
  !------------------------------------------!
  do c = 1, phi % pnt_grid % n_cells

    ! Under-relax central coefficient
    A % val(A % dia(c)) = A % val(A % dia(c)) / phi % urf

    ! Under-relax central coefficient
    b(c) = b(c) + A % val(A % dia(c)) * (1.0 - phi % urf) * phi % n(c)

  end do

  end subroutine
