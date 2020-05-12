!==============================================================================!
  subroutine Cgns_Mod_Initialize_Counters
!------------------------------------------------------------------------------!
!   Initializes counters used while CGNS file is being read.                   !
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  cnt_cells    = 0  ! number of cells (except boundary and interface)

  cnt_hex = 0  ! number of hexahedral  cells
  cnt_pyr = 0  ! number of pyramid     cells
  cnt_wed = 0  ! number of wedge/prism cells
  cnt_tet = 0  ! number of tetrahedral cells

  end subroutine
