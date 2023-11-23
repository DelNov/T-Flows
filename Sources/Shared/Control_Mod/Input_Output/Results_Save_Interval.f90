!==============================================================================!
  subroutine Results_Save_Interval(Control, val, verbose)
!------------------------------------------------------------------------------!
!>  Reads results save interval from the control file.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control  !! parent class
  integer             :: val      !! results save interval
  logical,   optional :: verbose  !! controls output verbosity
!==============================================================================!

  call Control % Read_Int_Item('RESULTS_SAVE_INTERVAL', 12, val, verbose)

  end subroutine
