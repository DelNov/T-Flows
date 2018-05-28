!==============================================================================!
  subroutine Cgns_Mod_Initialize_Counters
!------------------------------------------------------------------------------!
!   Initiallizes counters used while CGNS file is being read.                  !
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  cnt_nodes     = 0
  cnt_cells     = 0
  cnt_blocks    = 0
  cnt_bnd_conds = 0
  cnt_bnd_conds = 0

  cnt_hex = 0
  cnt_pyr = 0
  cnt_wed = 0
  cnt_tet = 0
  cnt_tri = 0
  cnt_qua = 0

  cnt_x = 0
  cnt_y = 0
  cnt_z = 0

  end subroutine
