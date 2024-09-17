!==============================================================================!
  subroutine Sca_O_Dia(Lin, n, b, s, Acon, Aval)
!------------------------------------------------------------------------------!
!>  Front-end for calculation scalar over matrix diagonal operation.
!------------------------------------------------------------------------------!
!   Note: Using intent clause here, was causing slower runs and crashes        !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Linalg_Type)            :: Lin   !! parent class
  integer, intent(in)           :: n     !! size of vectors
  real                          :: b(n)  !! result vector
  real                          :: s     !! scalar
  type(Sparse_Con_Type), target :: Acon  !! operand connectivity matrix
  type(Sparse_Val_Type), target :: Aval  !! operand values matrix
!-----------------------------------[Locals]-----------------------------------!
  integer, pointer :: a_dia(:)
  real,    pointer :: a_val(:)
  integer          :: i, nz
!==============================================================================!

  ! Take aliases
  nz    =  Acon % nonzeros
  a_dia => Acon % dia
  a_val => Aval % val

  !$acc parallel loop independent &
  !$acc present(  &
  !$acc   b,  &
  !$acc   a_val,  &
  !$acc   a_dia   &
  !$acc )
  do i = 1, n
    b(i) = s / a_val(a_dia(i))
  end do
  !$acc end parallel

  end subroutine

