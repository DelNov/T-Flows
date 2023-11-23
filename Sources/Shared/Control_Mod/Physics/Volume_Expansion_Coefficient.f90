!==============================================================================!
  subroutine Volume_Expansion_Coefficient(Control, cor, verbose)
!------------------------------------------------------------------------------!
!>  Reads the value of volume expansion coefficient from the control file.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control  !! parent class
  real,   intent(out) :: cor      !! volume expansion coefficient
  logical,   optional :: verbose  !! controls output verbosity
!-----------------------------------[Locals]-----------------------------------!
  real :: def
!==============================================================================!

  def = 1.0

  call Control % Read_Real_Item('VOLUME_EXPANSION_COEFFICIENT',  &
                                 def, cor, verbose)

  end subroutine
