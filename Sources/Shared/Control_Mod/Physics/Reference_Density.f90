!==============================================================================!
  subroutine Reference_Density(Control, d_ref, verbose)
!------------------------------------------------------------------------------!
!>  Reads the value of reference density from the control file.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control  !! parent class
  real,   intent(out) :: d_ref    !! value of reference density
  logical,   optional :: verbose  !! controls output verbosity
!-----------------------------------[Locals]-----------------------------------!
  real :: def
!==============================================================================!

  def = 0.0

  call Control % Read_Real_Item('REFERENCE_DENSITY', def, d_ref, verbose)

  end subroutine
