!==============================================================================!
  subroutine Numerics_Mod_Under_Relax(phi, a, b)
!------------------------------------------------------------------------------!
!   Purpose: Under-relax system of equations before calling linear solver.     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Var_Type)    :: phi
  type(Matrix_Type) :: a
  real              :: b(:)
!-----------------------------------[Locals]-----------------------------------!
  integer :: c
!==============================================================================!

  ! Browse through cells and under-relax
  do c = 1, phi % pnt_grid % n_cells
    b(c) = b(c)                                                  &
         + a % val(a % dia(c)) * (1.0 - phi % urf) * phi % n(c)  &
         / phi % urf
    a % val(a % dia(c)) = a % val(a % dia(c))  &
                        / phi % urf
  end do

  end subroutine
