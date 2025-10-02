!==============================================================================!
  subroutine Latent_Heat(Control, val, verbose)
!------------------------------------------------------------------------------!
!>  Reads the value of latent from the control file.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control  !! parent class
  real,   intent(out) :: val      !! latent heat value
  logical,   optional :: verbose  !! controls output verbosity
!==============================================================================!

  call Control % Read_Real_Item('LATENT_HEAT', 1.0, val, verbose)

  end subroutine
