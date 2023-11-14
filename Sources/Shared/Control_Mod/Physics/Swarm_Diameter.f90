!==============================================================================!
  subroutine Swarm_Diameter(Control, s_dia, verbose)
!------------------------------------------------------------------------------!
!>  Reads swarm diameter from the control file.  (This assumes that all
!>  particles have the same size.  For non-uniform distributions resort
!>  to user functions.)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control  !! parent class
  real,   intent(out) :: s_dia    !! particles' diameters
  logical,   optional :: verbose  !! controls output verbosity
!-----------------------------------[Locals]-----------------------------------!
  real :: def
!==============================================================================!

  def = 2.5e-5

  call Control % Read_Real_Item('SWARM_DIAMETER', def, s_dia, verbose)

  end subroutine
