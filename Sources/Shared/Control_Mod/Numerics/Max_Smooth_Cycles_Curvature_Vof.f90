!==============================================================================!
  subroutine Max_Smooth_Cycles_Curvature_Vof(Control, val, verbose)
!------------------------------------------------------------------------------!
!>  Reads maximum number of cycles for curvature estimation in VOF.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type)  :: Control  !! parent class
  integer, intent(out) :: val      !! max smoothing cycles for curvature in VOF
  logical, optional    :: verbose  !! controls output verbosity
!==============================================================================!

  call Control % Read_Int_Item('MAX_SMOOTHING_CYCLES_CURVATURE_VOF',  &
                                2, val, verbose)

  end subroutine
