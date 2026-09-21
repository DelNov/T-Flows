!==============================================================================!
  subroutine Buoyancy_Cutoff_Height(Control, cutoff_z, verbose)
!------------------------------------------------------------------------------!
!>  Reads, from the control file, an optional height above which the
!>  buoyancy force is set to zero.  Meant as a "sponge" for a
!>  computational domain that extends above the physically modelled
!>  region (e.g. a lab tank truncated for meshing convenience).  If not
!>  set, defaults to HUGE, so the cutoff never triggers.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control    !! parent class
  real,   intent(out) :: cutoff_z  !! height above which buoyancy is zeroed
  logical,   optional :: verbose    !! controls output verbosity
!-----------------------------------[Locals]-----------------------------------!
  real :: def
!==============================================================================!

  def = HUGE

  call Control % Read_Real_Item('BUOYANCY_CUTOFF_HEIGHT', def, cutoff_z,  &
                                verbose)

  end subroutine
