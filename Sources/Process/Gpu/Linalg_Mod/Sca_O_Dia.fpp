!==============================================================================!
  subroutine Sca_O_Dia(Lin, n, b, s, A)
!------------------------------------------------------------------------------!
!>  Front-end for calculation scalar over matrix diagonal operation.
!------------------------------------------------------------------------------!
!   Note: Using intent clause here, was causing slower runs and crashes        !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Linalg_Type)        :: Lin   !! parent class
  integer, intent(in)       :: n     !! size of vectors
  real                      :: b(n)  !! result vector
  real                      :: s     !! scalar
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
    b(i) = s / a_val(a_dia(i))
  end do
  !$tf-acc loop end

  end subroutine

