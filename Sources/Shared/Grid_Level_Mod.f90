!==============================================================================!
  module Grid_Level_Mod
!------------------------------------------------------------------------------!
!   This module is to facilitate the algebraic multigrid method in the code.   !
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !---------------------!
  !                     !
  !   Grid level type   !
  !                     !
  !---------------------!
  type Grid_Level_Type

    ! Number of cells and faces at each level
    integer :: n_cells
    integer :: n_faces

    ! Cell and face numbers from current to coarser levels
    integer, allocatable :: cell(:)
    integer, allocatable :: face(:)

    integer, allocatable :: n_finest_cells(:)

    ! Cell and face at coarser level
!   integer, allocatable :: cell_at_coarser(:)
!   integer, allocatable :: face_at_coarser(:)

    ! Faces' neigboring (surrounding) cells
    integer, allocatable :: faces_c(:,:)

    ! Cell coordinates
    real, allocatable :: xc(:), yc(:), zc(:)
  end type

! contains

! include 'Grid_Mod/Allocate_Levels.f90'

  end module
