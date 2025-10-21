!==============================================================================!
  real function Normalized_Root_Mean_Square(Lin, n, b, A)
!------------------------------------------------------------------------------!
!>  Front-end for calculation sqrt(vector^2 over matrix diagonal^2 operation).
!------------------------------------------------------------------------------!
!   Note: Using intent clause here, was causing slower runs and crashes        !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Linalg_Type)        :: Lin   !! parent class
  integer, intent(in)       :: n     !! size of vectors
  real                      :: b(n)  !! vector
  type(Sparse_Type), target :: A     !! matrix
!-----------------------------------[Locals]-----------------------------------!
  integer, pointer :: a_dia(:)
  real,    pointer :: a_val(:)
  real             :: rms
  integer          :: i, nz
!==============================================================================!

  ! Take aliases
  nz    =  A % nonzeros
  a_dia => A % dia
  a_val => A % val

  rms = 0.0

  !$tf-acc loop begin
  do i = 1, n
    rms = rms + b(i)**2 / a_val(a_dia(i))**2
  end do
  !$tf-acc loop end

  call Global % Sum_Real(rms)

  Normalized_Root_Mean_Square = sqrt(rms)

  end function

