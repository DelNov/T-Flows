!==============================================================================!
  subroutine Gradient_Method_For_Wall_Distance(Control, scheme_name, verbose)
!------------------------------------------------------------------------------!
!   Reading gradient method for wall distance                                  !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type)        :: Control
  character(SL), intent(out) :: scheme_name
  logical, optional          :: verbose
!==============================================================================!

  call Control % Read_Char_Item('GRADIENT_METHOD_FOR_WALL_DISTANCE',  &
                                'gauss_theorem',                      &
                                 scheme_name, verbose)
  call String % To_Upper_Case(scheme_name)

  end subroutine
