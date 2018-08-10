!==============================================================================!
  subroutine Cgns_Mod_Initialize_Counters
!------------------------------------------------------------------------------!
!   Initializes counters used while CGNS file is being read.                   !
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  cnt_nodes     = 0 ! number of nodes
  cnt_cells     = 0 ! number of cells (except boundary and interface)
  cnt_blocks    = 0 ! number of block
  cnt_bnd_conds = 0 ! number of boundary

  cnt_hex = 0 ! number of hexahedral  cells
  cnt_pyr = 0 ! number of pyramid     cells
  cnt_wed = 0 ! number of wedge/prism cells
  cnt_tet = 0 ! number of tetrahedral cells

  cnt_bnd_tri = 0 ! number of triangles defining boundary conditions
  cnt_bnd_qua = 0 ! number of quads     defining boundary conditions

  cnt_int     = 0 ! number of unique interfaces
  cnt_int_tri = 0 ! number of triangles defining block interface
  cnt_int_qua = 0 ! number of quads     defining block interface

  end subroutine