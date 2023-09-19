!==============================================================================!
  subroutine Create_Native(Nat, Grid)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Native_Type),       intent(out) :: Nat
  type(Grid_Type),  target, intent(in)  :: Grid
!==============================================================================!

  Nat % pnt_grid => Grid
  Assert(associated(Nat % pnt_grid))

  if(First_Proc()) print *, '# Determining matrix topology.'

  call Nat % A % Create_Matrix(Grid)
  call Nat % M % Create_Matrix(Grid)
  call Vector_Mod_Allocate(Nat % b, Grid)

  if(First_Proc()) print *, '# Finished !'

  end subroutine
