!==============================================================================!
  subroutine Control_Mod_Swarm_Subgrid_Scale_Model(val, verbose)
!------------------------------------------------------------------------------!
!   Reading swarm SGS model from the control file.                             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  character(SL), intent(out) :: val
  logical,       optional    :: verbose
!==============================================================================!

  call Control_Mod_Read_Char_Item('SWARM_SUBGRID_SCALE_MODEL', & 
                                  'none', val, verbose)

  end subroutine
