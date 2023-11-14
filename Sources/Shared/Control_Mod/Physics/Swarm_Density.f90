!==============================================================================!
  subroutine Swarm_Density(Control, s_den, verbose)
!------------------------------------------------------------------------------!
!>  Reads swarm density from the control file.  (This assumes that all
!>  particles have the same density, which will probably hold in most cases.)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control  !! parent class
  real,   intent(out) :: s_den    !! swarm (particles') mass density
  logical,   optional :: verbose  !! controls output verbosity
!-----------------------------------[Locals]-----------------------------------!
  real :: def
!==============================================================================!

  def = 1000.0

  call Control % Read_Real_Item('SWARM_DENSITY', def, s_den, verbose)

  end subroutine
