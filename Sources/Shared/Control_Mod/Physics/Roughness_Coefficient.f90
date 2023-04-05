!==============================================================================!
  subroutine Roughness_Coefficient(Control, val, verbose)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control
  real,   intent(out) :: val(:)
  logical,   optional :: verbose
!-----------------------------------[Locals]-----------------------------------!
  real :: con
!==============================================================================!

  call Control % Read_Real_Item('ROUGHNESS_COEFFICIENT', 0.0, con, verbose)

  ! Set the same value everywhere
  val(:) = con

  end subroutine
