!==============================================================================!
  subroutine Sys_Normalize(Lin, n, fn, Acon, Aval, b)
!------------------------------------------------------------------------------!
!>  Front-end for scaling (normalizing) a linear system of equations
!------------------------------------------------------------------------------!
!   Note: Using intent clause here, was causing slower runs and crashes        !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Linalg_Type)    :: Lin   !! parent class
  integer,  intent(in)  :: n     !! size of vectors
  real,     intent(out) :: fn    !! factor of normalization
  type(Sparse_Con_Type) :: Acon  !! operand connectivity matrix
  type(Sparse_Val_Type) :: Aval  !! operand values matrix
  real                  :: b(n)  !! right hand side vector
!-----------------------------------[Locals]-----------------------------------!
  integer :: nz
!==============================================================================!

  ! Take aliases
  nz = Acon % nonzeros

  ! ... and then compute matrix vector product
  ! on the device attached to this processor
  call Lin % Sys_Normalize_Acc(n,             &
                               nz,            &
                               fn,            &
                               Aval % val,    &
                               Acon % dia,    &
                               Aval % d_inv,  &
                               b)

  end subroutine

