!==============================================================================!
  subroutine Swarm_Subgrid_Scale_Model(Control, val, verbose)
!------------------------------------------------------------------------------!
!   Reading swarm SGS model from the control file.                             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type)        :: Control
  character(SL), intent(out) :: val
  logical,       optional    :: verbose
!==============================================================================!

  call Control % Read_Char_Item('SWARM_SUBGRID_SCALE_MODEL', &
                                'none', val, verbose)

  end subroutine
