!==============================================================================!
  subroutine Swarm_Coefficient_Of_Restitution(Control, cor, verbose)
!------------------------------------------------------------------------------!
!>  Reads swarm coefficient of restitution from the control file.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control  !! parent class
  real,   intent(out) :: cor      !! value of the coefficient of restitution
  logical,   optional :: verbose  !! controls output verbosity
!-----------------------------------[Locals]-----------------------------------!
  real :: def
!==============================================================================!

  def = 1.0

  call Control % Read_Real_Item('SWARM_COEFFICIENT_OF_RESTITUTION',  &
                                 def, cor, verbose)

  end subroutine
