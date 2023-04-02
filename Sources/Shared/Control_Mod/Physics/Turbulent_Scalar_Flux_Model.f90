!==============================================================================!
  subroutine Control_Mod_Turbulent_Scalar_Flux_Model(val, verbose)
!------------------------------------------------------------------------------!
!   Reading turbulent heat flux model from the control file.                   !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  character(SL), intent(out) :: val
  logical,       optional    :: verbose
!==============================================================================!

  call Control % Read_Char_Item('TURBULENT_SCALAR_FLUX_MODEL', 'SGDH',  &
                                 val, verbose)
  call String % To_Upper_Case(val)

  end subroutine
