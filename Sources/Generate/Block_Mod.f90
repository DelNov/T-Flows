!==============================================================================!
  module Block_Mod
!------------------------------------------------------------------------------!
!>  This module is responsible for representing the concept of a hexahedral
!>  block within the context of mesh generation in the sub-program Generate.
!>  Together with Line_Mod (and its Line_Type) it helps do define a
!>  computational domain in Domain_Mod.
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !----------------!
  !   Block type   !
  !----------------!
  !> Encapsulates the properties and characteristics
  !> of a structured block within the mesh.
  type Block_Type

    integer :: n_nodes            !! number of nodes in the block
    integer :: n_cells            !! number of cells in the block
    integer :: faces(6,4)         !! six faces of a hexahedron and their nodes
    integer :: corners(0:8)       !! 0 is for orientation, 1-8 for points
    integer :: resolutions(3)     !! block resolution (ni, nj, nk)
    real    :: weights(3)         !! block weights for node clustering
    real    :: face_weights(6,3)  !! weights for each of the block's faces
                                  !! used for node clustering
  end type

  end module
