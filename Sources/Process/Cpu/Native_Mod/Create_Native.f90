!==============================================================================!
  subroutine Create_Native(Nat, Grid)
!------------------------------------------------------------------------------!
!>  Initializes a Native_Type object for linear solver operations. The
!>  subroutine associates the Native object with a given numerical grid and
!>  prepares matrices and vectors for linear system computations.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Native_Type),       intent(out) :: Nat   !! parent, Native_Type object
  type(Grid_Type),  target, intent(in)  :: Grid  !! grid on which it's defined
!==============================================================================!

  Nat % pnt_grid => Grid
  Assert(associated(Nat % pnt_grid))

  if(First_Proc()) print *, '# Determining matrix topology.'

  call Nat % A % Create_Matrix(Grid)
  call Nat % M % Create_Matrix(Grid)
  call Vector_Mod_Allocate(Nat % b, Grid)

  if(First_Proc()) print *, '# Finished !'

  end subroutine
