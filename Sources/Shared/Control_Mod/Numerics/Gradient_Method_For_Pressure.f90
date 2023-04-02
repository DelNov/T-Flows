!==============================================================================!
  subroutine Control_Mod_Gradient_Method_For_Pressure(scheme_name, verbose)
!------------------------------------------------------------------------------!
!   Reading gradient method for pressure.                                      !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  character(SL), intent(out) :: scheme_name
  logical, optional          :: verbose
!==============================================================================!

  call Control % Read_Char_Item('GRADIENT_METHOD_FOR_PRESSURE',  &
                                'least_squares',                 &
                                 scheme_name, verbose)
  call String % To_Upper_Case(scheme_name)

  end subroutine
