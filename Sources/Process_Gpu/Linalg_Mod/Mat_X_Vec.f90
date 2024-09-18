!==============================================================================!
  subroutine Mat_X_Vec(Lin, n, c, Acon, Aval, b)
!------------------------------------------------------------------------------!
!>  Front-end for calculation of sparse-matrix vector multiplication.
!------------------------------------------------------------------------------!
!   Note: Using intent clause here, was causing slower runs and crashes        !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Linalg_Type)            :: Lin   !! parent class
  integer, intent(in)           :: n     !! size of vectors
  real                          :: c(n)  !! result vector
  type(Sparse_Con_Type), target :: Acon  !! operand connectivity matrix
  type(Sparse_Val_Type), target :: Aval  !! operand values matrix
  real                          :: b(n)  !! operand vector
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
  real,    pointer         :: a_val(:)
  integer, pointer         :: a_col(:), a_row(:)
  real                     :: temp
  integer                  :: nz, i, j, ij
!==============================================================================!

  ! Take aliases
  Grid  => Acon % pnt_grid
  nz    =  Acon % nonzeros
  a_col => Acon % col
  a_row => Acon % row
  a_val => Aval % val

  ! Refresh the operand vector over processor buffers ...
  call Grid % Exchange_Inside_Cells_Real(b(1:n))

  !$acc parallel loop independent  &
  !$acc present(  &
  !$acc   a_row,  &
  !$acc   a_col,  &
  !$acc   c,  &
  !$acc   a_val,  &
  !$acc   b   &
  !$acc )
  do i = 1, n
    temp = 0.0

  !$acc loop seq
    do ij = a_row(i), a_row(i+1) - 1
      j = a_col(ij)
      temp = temp + a_val(ij) * b(j)
    end do
  !$acc end loop

    c(i) = temp
  end do
  !$acc end parallel

  end subroutine

