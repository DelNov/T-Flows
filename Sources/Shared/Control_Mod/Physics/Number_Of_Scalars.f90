!==============================================================================!
  subroutine Number_Of_Scalars(Control, val, verbose)
!------------------------------------------------------------------------------!
!>  Reads the number of simulated scalar transport equations.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type)  :: Control  !! parent class
  integer, intent(out) :: val      !! number of scalars
  logical, optional    :: verbose  !! controls output verbosity
!==============================================================================!

  call Control % Read_Int_Item('NUMBER_OF_SCALARS', 0, val, verbose)

  end subroutine
