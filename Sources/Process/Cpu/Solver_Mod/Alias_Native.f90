!==============================================================================!
  subroutine Alias_Native(Sol, A, b)
!------------------------------------------------------------------------------!
!>  This subroutine is used for streamlining the access to the system matrix
!>  and the right-hand side vector when using native solvers. Its primary
!>  purpose is to create aliases for these components.  In practice, once
!>  Alias_Native is called, subsequent code sections can refer to the system
!>  matrix and right-hand side vector using the simple and short 'A' and 'b'
!>  notation instead of the longer form: Sol % Nat % A and Sol % Nat % b.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Solver_Type), target  :: Sol   !! parent class
  type(Matrix_Type),  pointer :: A     !! alias (abbreviation) to Sol % Nat % A
  real,               pointer :: b(:)  !! alias (abbreviation) to Sol % Nat % b
!==============================================================================!

  A => Sol % Nat % A
  b => Sol % Nat % b % val

  end subroutine
