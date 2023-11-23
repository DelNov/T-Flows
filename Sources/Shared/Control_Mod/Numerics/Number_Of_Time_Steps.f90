!==============================================================================!
  subroutine Number_Of_Time_Steps(Control, val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type)  :: Control  !! parent class
  integer, intent(out) :: val      !! number of time steps
  logical, optional    :: verbose  !! controls output verbosity
!==============================================================================!

  call Control % Read_Int_Item('NUMBER_OF_TIME_STEPS', 1200, val, verbose)

  end subroutine
