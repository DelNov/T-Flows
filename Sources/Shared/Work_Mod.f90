!==============================================================================!
  module Work_Mod
!------------------------------------------------------------------------------!
!   This is used only to store working variables to used instead of temporary  !
!   arrays throughout the code.  By doing so, memory is saved (reused),        !
!   as well as time for allocation and de-allocation.                          !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !-------------------------------!
  !   Cell-based working arrays   !
  !-------------------------------!

  ! Real; 01 - 10
  real, allocatable :: r_cell_01(:), r_cell_02(:), r_cell_03(:), r_cell_04(:)
  real, allocatable :: r_cell_05(:), r_cell_06(:), r_cell_07(:), r_cell_08(:)
  real, allocatable :: r_cell_09(:), r_cell_10(:)
  ! Real; 11 - 20
  real, allocatable :: r_cell_11(:), r_cell_12(:), r_cell_13(:), r_cell_14(:)
  real, allocatable :: r_cell_15(:), r_cell_16(:), r_cell_17(:), r_cell_18(:)
  real, allocatable :: r_cell_19(:), r_cell_20(:)
  ! Real; 21 - 30
  real, allocatable :: r_cell_21(:), r_cell_22(:), r_cell_23(:), r_cell_24(:)
  real, allocatable :: r_cell_25(:), r_cell_26(:), r_cell_27(:), r_cell_28(:)
  real, allocatable :: r_cell_29(:), r_cell_30(:)

  ! Integer; 01 - 10
  integer, allocatable :: i_cell_01(:), i_cell_02(:), i_cell_03(:)
  integer, allocatable :: i_cell_04(:), i_cell_05(:), i_cell_06(:)
  integer, allocatable :: i_cell_07(:), i_cell_08(:), i_cell_09(:)
  integer, allocatable :: i_cell_10(:)
  ! Integer; 11 - 20
  integer, allocatable :: i_cell_11(:), i_cell_12(:), i_cell_13(:)
  integer, allocatable :: i_cell_14(:), i_cell_15(:), i_cell_16(:)
  integer, allocatable :: i_cell_17(:), i_cell_18(:), i_cell_19(:)
  integer, allocatable :: i_cell_20(:)
  ! Integer; 21 - 30
  integer, allocatable :: i_cell_21(:), i_cell_22(:), i_cell_23(:)
  integer, allocatable :: i_cell_24(:), i_cell_25(:), i_cell_26(:)
  integer, allocatable :: i_cell_27(:), i_cell_28(:), i_cell_29(:)
  integer, allocatable :: i_cell_30(:)

  !-------------------------------!
  !   Face-based working arrays   !
  !-------------------------------!

  ! Real; 01 - 10
  real, allocatable :: r_face_01(:), r_face_02(:), r_face_03(:), r_face_04(:)
  real, allocatable :: r_face_05(:), r_face_06(:), r_face_07(:), r_face_08(:)
  real, allocatable :: r_face_09(:), r_face_10(:)
  ! Real; 11 - 20
  real, allocatable :: r_face_11(:), r_face_12(:), r_face_13(:), r_face_14(:)
  real, allocatable :: r_face_15(:), r_face_16(:), r_face_17(:), r_face_18(:)
  real, allocatable :: r_face_19(:), r_face_20(:)
  ! Real; 21 - 30
  real, allocatable :: r_face_21(:), r_face_22(:), r_face_23(:), r_face_24(:)
  real, allocatable :: r_face_25(:), r_face_26(:), r_face_27(:), r_face_28(:)
  real, allocatable :: r_face_29(:), r_face_30(:)

  ! Integer; 01 - 10
  integer, allocatable :: i_face_01(:), i_face_02(:), i_face_03(:)
  integer, allocatable :: i_face_04(:), i_face_05(:), i_face_06(:)
  integer, allocatable :: i_face_07(:), i_face_08(:), i_face_09(:)
  integer, allocatable :: i_face_10(:)
  ! Integer; 11 - 20
  integer, allocatable :: i_face_11(:), i_face_12(:), i_face_13(:)
  integer, allocatable :: i_face_14(:), i_face_15(:), i_face_16(:)
  integer, allocatable :: i_face_17(:), i_face_18(:), i_face_19(:)
  integer, allocatable :: i_face_20(:)
  ! Integer; 21 - 30
  integer, allocatable :: i_face_21(:), i_face_22(:), i_face_23(:)
  integer, allocatable :: i_face_24(:), i_face_25(:), i_face_26(:)
  integer, allocatable :: i_face_27(:), i_face_28(:), i_face_29(:)
  integer, allocatable :: i_face_30(:)

  ! Real; 01 - 08
  real, allocatable :: r_node_01(:), r_node_02(:), r_node_03(:), r_node_04(:)
  real, allocatable :: r_node_05(:), r_node_06(:), r_node_07(:), r_node_08(:)
  ! Integer; 01 - 08
  integer, allocatable :: i_node_01(:), i_node_02(:), i_node_03(:), i_node_04(:)
  integer, allocatable :: i_node_05(:), i_node_06(:), i_node_07(:), i_node_08(:)

  contains

  include 'Work_Mod/Allocate.f90'
  include 'Work_Mod/Allocate_Integer_Cells.f90'
  include 'Work_Mod/Allocate_Integer_Faces.f90'
  include 'Work_Mod/Allocate_Integer_Nodes.f90'
  include 'Work_Mod/Allocate_Real_Cells.f90'
  include 'Work_Mod/Allocate_Real_Faces.f90'
  include 'Work_Mod/Allocate_Real_Nodes.f90'

  end module
