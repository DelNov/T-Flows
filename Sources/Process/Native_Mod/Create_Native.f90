!==============================================================================!
  subroutine Create_Native(Nat, Grid)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Native_Type)       :: Nat
  type(Grid_Type),  target :: Grid
!==============================================================================!

  Nat % pnt_grid => Grid

  if(this_proc < 2) print *, '# Determining matrix topology.'

  call Nat % A % Create_Matrix(Grid)
  call Nat % M % Create_Matrix(Grid)
  call Nat % D % Create_Matrix(Grid)
  call Vector_Mod_Allocate(Nat % b, Grid)

  if(this_proc < 2) print *, '# Finished !'

  end subroutine
