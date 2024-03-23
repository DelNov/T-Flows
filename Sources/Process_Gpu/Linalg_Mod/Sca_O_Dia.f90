!==============================================================================!
  subroutine Sca_O_Dia(Lin, n, b, s, Acon, Aval)
!------------------------------------------------------------------------------!
!>  Front-end for calculation scalar over matrix diagonal operation.
!------------------------------------------------------------------------------!
!   Note: Using intent clause here, was causing slower runs and crashes        !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Linalg_Type)    :: Lin   !! parent class
  integer, intent(in)   :: n     !! size of vectors
  real                  :: b(n)  !! result vector
  real                  :: s     !! scalar
  type(Sparse_Con_Type) :: Acon  !! operand connectivity matrix
  type(Sparse_Val_Type) :: Aval  !! operand values matrix
!-----------------------------------[Locals]-----------------------------------!
  integer :: nz
!==============================================================================!

  ! Take aliases
  nz = Acon % nonzeros

  ! ... and then compute matrix vector product
  ! on the device attached to this processor
  call Lin % Sca_O_Dia_Acc(n,           &
                           nz,          &
                           b,           &
                           s,           &
                           Aval % val,  &
                           Acon % dia)

  end subroutine

