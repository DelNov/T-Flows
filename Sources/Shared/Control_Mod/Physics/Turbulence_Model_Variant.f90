!==============================================================================!
  subroutine Turbulence_Model_Variant(Control, val, verbose)
!------------------------------------------------------------------------------!
!>  Reads turbulence model variant from the control file.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type)        :: Control  !! parent class
  character(SL), intent(out) :: val      !! variant of the turbulence model
  logical,       optional    :: verbose  !! controls output verbosity
!==============================================================================!

  call Control % Read_Char_Item('TURBULENCE_MODEL_VARIANT', 'stabilized',  &
                                 val, verbose)
  call String % To_Upper_Case(val)

  end subroutine
