!==============================================================================!
  subroutine Max_Substep_Cycles_Vof(Control, val, verbose)
!------------------------------------------------------------------------------!
!>  Reads maximum number of sub-step cycles in VOF.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type)  :: Control  !! parent class
  integer, intent(out) :: val      !! max sub-step cycles in VOF
  logical, optional    :: verbose  !! controls output verbosity
!==============================================================================!

  call Control % Read_Int_Item('MAX_SUBSTEP_CYCLES_VOF',  &
                                100, val, verbose)

  end subroutine
