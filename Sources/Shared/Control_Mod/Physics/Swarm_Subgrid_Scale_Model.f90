!==============================================================================!
  subroutine Control_Mod_Swarm_Subgrid_Scale_Model(s_sgs, verbose)
!------------------------------------------------------------------------------!
!             Reading swarm SGS model from the control file.                   !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  character(len=80) :: s_sgs
  logical, optional :: verbose
!==============================================================================!

  call Control_Mod_Read_Char_Item('SWARM_SUBGRID_SCALE_MODEL', & 
                                  'BROWNIAN_FUKAGATA', s_sgs, verbose)

  end subroutine
