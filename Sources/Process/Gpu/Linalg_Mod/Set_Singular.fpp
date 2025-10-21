!==============================================================================!
  subroutine Set_Singular(Lin, n, A)
!------------------------------------------------------------------------------!
!>  The Set_Singular subroutine in the Solver_Mod module plays an important
!>  role in preparing the solver for handling singular matrices. This
!>  functionality is especially relevant when dealing with matrices that might
!>  be deficient in rank or lack full invertibility, which is typically the
!>  case for the pressure correction solution.
!------------------------------------------------------------------------------!
!   Note: Using intent clause here, was causing slower runs and crashes        !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Linalg_Type)        :: Lin   !! parent class
  integer, intent(in)       :: n     !! size of vectors
  type(Sparse_Type), target :: A     !! operand matrix
!-----------------------------------[Locals]-----------------------------------!
  integer, pointer :: a_dia(:)
  real,    pointer :: a_val(:)
  integer          :: i, nz
!==============================================================================!

  ! Take aliases
  nz    =  A % nonzeros
  a_dia => A % dia
  a_val => A % val

  !$tf-acc loop begin
  do i = 1, n
    a_val(a_dia(i)) = a_val(a_dia(i)) * (1.0 + MICRO)
  end do
  !$tf-acc loop end

  end subroutine

