!==============================================================================!
  module Gen_Mod
!------------------------------------------------------------------------------!
!   Global variable definitions for the mesh generator.                        !
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  ! Twin nodes
  integer, allocatable :: twin_n(:,:)

  ! For each face, which neighbours are to each other cells which meet there
  ! The name derives from "face cell to cell"
  integer, allocatable :: face_c_to_c(:,:)

  ! Periodic and copy boundaries
  integer              :: n_periodic_cond
  integer, allocatable ::   periodic_cond(:,:)

  contains

#   include "Gen_Mod/Are_Nodes_Twins.f90"

  end module
