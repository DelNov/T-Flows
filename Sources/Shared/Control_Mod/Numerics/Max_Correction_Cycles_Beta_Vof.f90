!==============================================================================!
  subroutine Max_Correction_Cycles_Beta_Vof(Control, val, verbose)
!------------------------------------------------------------------------------!
!>  Reads maximum number of correction cycles for beta (in VOF).
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type)  :: Control  !! parent class
  integer, intent(out) :: val      !! correction cycles for beta VOF
  logical, optional    :: verbose  !! controls output verbosity
!==============================================================================!

  call Control % Read_Int_Item('MAX_CORRECTION_CYCLES_BETA_VOF',  &
                                2, val, verbose)

  end subroutine
