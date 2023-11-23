!==============================================================================!
  subroutine Gradient_Method_For_Pressure(Control, grad_method, verbose)
!------------------------------------------------------------------------------!
!>  Reads gradient method for pressure from control file.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type)        :: Control      !! parent class
  character(SL), intent(out) :: grad_method  !! gradient method
  logical, optional          :: verbose      !! controls output verbosity
!==============================================================================!

  call Control % Read_Char_Item('GRADIENT_METHOD_FOR_PRESSURE',  &
                                'least_squares',                 &
                                 grad_method, verbose)
  call String % To_Upper_Case(grad_method)

  end subroutine
