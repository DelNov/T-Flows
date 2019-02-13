!==============================================================================!
  subroutine Control_Mod_Turbulence_Model_Variant(val, verbose)
!------------------------------------------------------------------------------!
!   Reading turbulence model variant from the control file                     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  logical, optional :: verbose
!-----------------------------------[Locals]-----------------------------------!
  character(len=80) :: val
!==============================================================================!

  call Control_Mod_Read_Char_Item('TURBULENCE_MODEL_VARIANT', 'stabilized',  &
                                   val, verbose)
  call To_Upper_Case(val)

  end subroutine
