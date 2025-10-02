!==============================================================================!
  subroutine Prec_Form(Nat, Acon, Aval)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Native_Type),    target, intent(in) :: Nat   !! parent class
  type(Sparse_Con_Type), target, intent(in) :: Acon  !! connectivity matrix
  type(Sparse_Val_Type), target, intent(in) :: Aval  !! value matrix
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: Grid
  real, contiguous, pointer :: d_inv(:)
  integer                   :: i, ij
!==============================================================================!

  ! Take aliases
  Grid  => Acon % pnt_grid
  d_inv => Aval % d_inv

  ! Prepare matrix for diagonal preconditioning
  do i = 1, Grid % n_cells
    do ij = Acon % row(i), Acon % row(i+1) - 1
      if(ij .eq. Acon % dia(i)) then        ! store reciprocal of diagonal
        d_inv(i) = 1.0 / Aval % val(ij)
      end if
    end do
  end do

  end subroutine
