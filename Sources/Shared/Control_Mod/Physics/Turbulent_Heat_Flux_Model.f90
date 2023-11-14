!==============================================================================!
  subroutine Turbulent_Heat_Flux_Model(Control, val, verbose)
!------------------------------------------------------------------------------!
!>  Reads turbulence heat flux model (SGDH, GGDH) from the control file.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type)        :: Control  !! parent class
  character(SL), intent(out) :: val      !! turbulent heat flux model
  logical,       optional    :: verbose  !! controls output verbosity
!==============================================================================!

  call Control % Read_Char_Item('TURBULENT_HEAT_FLUX_MODEL', 'SGDH',  &
                                 val, verbose)
  call String % To_Upper_Case(val)

  end subroutine
