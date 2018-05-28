!==============================================================================!
  module Block_Mod
!------------------------------------------------------------------------------!
!   Hexahedral blocks making a domain used in "Generator"                      !
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !----------------!
  !   Block type   !
  !----------------!
  type Block_Type

    integer :: n_nodes            ! number of nodes in the block
    integer :: n_cells            ! number of cells in the block
    integer :: faces(6,4)         ! faces' nodes
    integer :: corners(0:8)       ! 0 is for orientation, 1-8 for points                
    integer :: resolutions(3)     ! ni, nj, nk
    real    :: weights(3)         ! block weights for node clustering
    real    :: face_weights(6,3)  ! block weights for node clustering

  end type

  end module
