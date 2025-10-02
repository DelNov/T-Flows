!==============================================================================!
  subroutine Number_Of_Phases(Control, val, verbose)
!------------------------------------------------------------------------------!
!>  Reads number of phases from the control file.  (This function might be
!>  obsolete, VOF simulations are limited to two phases.)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type)  :: Control  !! parent class
  integer, intent(out) :: val      !! number of phases
  logical, optional    :: verbose  !! controls output verbosity
!==============================================================================!

  call Control % Read_Int_Item('NUMBER_OF_PHASES', 1, val, verbose)

  end subroutine
