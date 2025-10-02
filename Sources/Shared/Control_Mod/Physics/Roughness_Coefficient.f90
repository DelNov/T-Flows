!==============================================================================!
  subroutine Roughness_Coefficient(Control, val, verbose)
!------------------------------------------------------------------------------!
!>  Reads wall roughness coefficient from the control file.  (This could be
!>  obsolete, see also Rough_Walls.)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control  !! parent class
  real,   intent(out) :: val(:)   !! roughness coefficient
  logical,   optional :: verbose  !! controls output verbosity
!-----------------------------------[Locals]-----------------------------------!
  real :: con
!==============================================================================!

  call Control % Read_Real_Item('ROUGHNESS_COEFFICIENT', 0.0, con, verbose)

  ! Set the same value everywhere
  val(:) = con

  end subroutine
