!==============================================================================!
  subroutine Prec_Form(Nat, A)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Native_Type), target, intent(in) :: Nat  !! parent, Native_Type object
  type(Sparse_Type),  target, intent(in) :: A    !! an existing matrix
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),   pointer :: Grid
  real, contiguous,  pointer :: d_inv(:)
  integer                    :: i, ij
!==============================================================================!

  ! Take aliases
  Grid  => A % pnt_grid
  d_inv => A % d_inv

  ! Prepare matrix for diagonal preconditioning
  do i = 1, Grid % n_cells
    do ij = A % row(i), A % row(i+1) - 1
      if(ij .eq. A % dia(i)) then        ! store reciprocal of diagonal
        d_inv(i) = 1.0 / A % val(ij)
      end if
    end do
  end do

  end subroutine
