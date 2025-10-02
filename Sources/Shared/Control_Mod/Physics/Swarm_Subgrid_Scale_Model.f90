!==============================================================================!
  subroutine Swarm_Subgrid_Scale_Model(Control, val, verbose)
!------------------------------------------------------------------------------!
!>  Reads swarm SGS model from the control file.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type)        :: Control  !! parent class
  character(SL), intent(out) :: val      !! SGS model for particles
  logical,       optional    :: verbose  !! controls output verbosity
!==============================================================================!

  call Control % Read_Char_Item('SWARM_SUBGRID_SCALE_MODEL', &
                                'none', val, verbose)

  end subroutine
