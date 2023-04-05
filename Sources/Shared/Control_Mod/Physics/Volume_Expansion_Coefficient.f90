!==============================================================================!
  subroutine Volume_Expansion_Coefficient(Control, cor, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control
  real,   intent(out) :: cor
  logical,   optional :: verbose
!-----------------------------------[Locals]-----------------------------------!
  real :: def
!==============================================================================!

  def = 1.0

  call Control % Read_Real_Item('VOLUME_EXPANSION_COEFFICIENT',  &
                                 def, cor, verbose)

  end subroutine
