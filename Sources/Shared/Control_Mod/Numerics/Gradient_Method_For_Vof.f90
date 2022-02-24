!==============================================================================!
  subroutine Control_Mod_Gradient_Method_For_Vof(scheme_name, verbose)
!------------------------------------------------------------------------------!
!   Reading gradient method for multiphase variables.                          !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  character(SL), intent(out) :: scheme_name
  logical, optional          :: verbose
!==============================================================================!

  call Control_Mod_Read_Char_Item('GRADIENT_METHOD_FOR_VOF',  &
                                  'least_squares',            &
                                   scheme_name, verbose)
  call To_Upper_Case(scheme_name)

  end subroutine
