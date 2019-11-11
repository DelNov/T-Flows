!==============================================================================!
  subroutine Control_Mod_Turbulence_Model(val, verbose)
!------------------------------------------------------------------------------!
!   Reading turbulence model from the control file.                            !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  character(len=80) :: val
  logical, optional :: verbose
!==============================================================================!

  call Control_Mod_Read_Char_Item('TURBULENCE_MODEL', 'none',  &
                                   val, verbose)
  call To_Upper_Case(val)

  end subroutine
