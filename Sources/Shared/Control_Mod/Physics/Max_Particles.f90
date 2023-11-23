!==============================================================================!
  subroutine Max_Particles(Control, val, verbose)
!------------------------------------------------------------------------------!
!>  Reads maximum number of particles from the control file.  (This might be
!>  obsolete/redundant, because it can be set through user functions.  Check!)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type)  :: Control  !! parent class
  integer, intent(out) :: val      !! maximum number of particles
  logical, optional    :: verbose  !! controls output verbosity
!==============================================================================!

  call Control % Read_Int_Item('MAX_PARTICLES', 0, val, verbose)

  end subroutine
