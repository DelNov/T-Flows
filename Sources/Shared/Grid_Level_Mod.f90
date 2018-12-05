!==============================================================================!
  module Grid_Level_Mod
!------------------------------------------------------------------------------!
!   This module is to facilitate the algebraic multigrid method in the code.   !
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !-----------------------------------------!
  !   Maximum number of multi-grid levels   !
  !-----------------------------------------!
  integer, parameter :: MAX_MG_LEVELS = 12

  !---------------------!
  !                     !
  !   Grid level type   !
  !                     !
  !---------------------!
  type Grid_Level_Type

    ! Number of cells and faces at each level
    integer :: n_cells
    integer :: n_faces

    ! Cell and face numbers from finest (1) to coarser levels (.gt. 1)
    integer, allocatable :: cell(:)
    integer, allocatable :: face(:)

    ! Cell at coarser level
    integer, allocatable :: coarser_c(:)

    integer, allocatable :: n_finest_cells(:)

    ! Faces' neigboring (surrounding) cells
    integer, allocatable :: faces_c(:,:)

    ! Cell coordinates
    real, allocatable :: xc(:), yc(:), zc(:)
  end type

! contains
!
! include 'Grid_Level_Mod/Allocate.f90'

  end module
