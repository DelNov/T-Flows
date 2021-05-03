!==============================================================================!
  subroutine Control_Mod_Gradient_Method_For_Energy(scheme_name, verbose)
!------------------------------------------------------------------------------!
!   Reading gradient method for energy (temperature, in practice)              !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  character(SL)     :: scheme_name
  logical, optional :: verbose
!==============================================================================!

  call Control_Mod_Read_Char_Item('GRADIENT_METHOD_FOR_ENERGY',  &
                                  'least_squares',               &
                                   scheme_name, verbose)
  call To_Upper_Case(scheme_name)

  end subroutine
