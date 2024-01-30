!==============================================================================!
  subroutine Lee_Model_Coefficients(Control, c_lee, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control
  real,   intent(out) :: c_lee(2)
  logical,   optional :: verbose
!-----------------------------------[Locals]-----------------------------------!
  real :: def(2)  ! defualt values
!==============================================================================!

  ! Set the default values (do these make sense?)
  data def / 0.0, 0.0 /

  ! Read from the control file
  call Control % Read_Real_Vector('LEE_MODEL_COEFFICIENTS',  &
                                  2, def, c_lee, verbose)

  end subroutine
