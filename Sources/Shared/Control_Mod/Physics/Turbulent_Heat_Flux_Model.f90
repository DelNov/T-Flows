!==============================================================================!
  subroutine Control_Mod_Turbulent_Heat_Flux_Model(val, verbose)
!------------------------------------------------------------------------------!
!   Reading turbulent heat flux model from the control file.                   !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  character(SL)     :: val
  logical, optional :: verbose
!==============================================================================!

  call Control_Mod_Read_Char_Item('TURBULENT_HEAT_FLUX_MODEL', 'SGDH',  &
                                   val, verbose)
  call To_Upper_Case(val)

  end subroutine
