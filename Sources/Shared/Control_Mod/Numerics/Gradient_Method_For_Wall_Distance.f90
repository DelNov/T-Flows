!==============================================================================!
  subroutine Control_Mod_Gradient_Method_For_Wall_Distance(scheme_name,  &
                                                           verbose)
!------------------------------------------------------------------------------!
!   Reading gradient method for wall distance                                  !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  character(SL), intent(out) :: scheme_name
  logical, optional          :: verbose
!==============================================================================!

  call Control % Read_Char_Item('GRADIENT_METHOD_FOR_WALL_DISTANCE',  &
                                'gauss_theorem',                      &
                                 scheme_name, verbose)
  call String % To_Upper_Case(scheme_name)

  end subroutine
