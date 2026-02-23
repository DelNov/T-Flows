!==============================================================================!
  subroutine Scalars_Diffusivities(Control, val, nsc, verbose)
!------------------------------------------------------------------------------!
!>  Reads the value of scalar diffusivity from the control file.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type)  :: Control           !! parent class
  real,    intent(out) :: val(MAX_SCALARS)  !! scalar diffusivity value
  integer,  intent(in) :: nsc               !! number of scalars
  logical,    optional :: verbose           !! controls output verbosity
!-----------------------------------[Locals]-----------------------------------!
  real :: def(MAX_SCALARS)
!==============================================================================!

  def = 1.0e-6

  call Control % Read_Real_Vector('SCALARS_DIFFUSIVITIES', nsc,  &
                                  def, val, verbose)

  end subroutine
