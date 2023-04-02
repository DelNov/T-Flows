!==============================================================================!
  subroutine Control_Mod_Turbulence_Model_Variant(val, verbose)
!------------------------------------------------------------------------------!
!   Reading turbulence model variant from the control file                     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  character(SL), intent(out) :: val
  logical,       optional    :: verbose
!==============================================================================!

  call Control % Read_Char_Item('TURBULENCE_MODEL_VARIANT', 'stabilized',  &
                                 val, verbose)
  call String % To_Upper_Case(val)

  end subroutine
