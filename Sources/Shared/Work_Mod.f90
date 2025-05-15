#include "../Shared/Assert.h90"

! Set this variable to 0, to avoid setting work array to zero every time you
! connect them.  This will not lead to any noticable differnce in performance,
! but not resetting connecting arrays to zero might help with user arrays, if
! a user doesn't want to forget the values in between calls to a function.
#define  RESET_VALUES_TO_ZERO  1

!==============================================================================!
  module Work_Mod
!------------------------------------------------------------------------------!
!>  The Work_Mod module in T-Flows software serves as a centralized resource
!>  manager for working variables, replacing the frequent use of temporary
!>  arrays throughout the code. This module aims to optimize memory usage by
!>  reusing memory spaces, thus saving both memory and time that would
!>  otherwise be spent in repeated allocation/deallocation of temporary arrays.
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  integer, parameter :: MAX_WORK_ARRAYS = 128

  !---------------------!
  !  Real arrays type   !
  !---------------------!
  !> Simple container of a real allocatable arrays.
  type Real_Arrays_Type
    real, allocatable :: array(:)
  end type

  !------------------------!
  !  Integer arrays type   !
  !------------------------!
  !> Simple container of a integer allocatable arrays.
  type Int_Arrays_Type
    integer, allocatable :: array(:)
  end type

  !--------------!
  !  Work type   !
  !--------------!
  !> The Work_Type within this module encapsulates various arrays for both
  !> integer and real data types, specifically tailored for different grid
  !> elements like cells, faces, and nodes.
  type Work_Type

    ! Real working arrays
    type(Real_Arrays_Type), allocatable :: r_cell(:)  !! cell-based work array
    type(Real_Arrays_Type), allocatable :: r_face(:)  !! face-based work array
    type(Real_Arrays_Type), allocatable :: r_node(:)  !! node-based work array

    ! Integer working arrays
    type(Int_Arrays_Type), allocatable :: i_cell(:)  !! cell-based work array
    type(Int_Arrays_Type), allocatable :: i_face(:)  !! face-based work array
    type(Int_Arrays_Type), allocatable :: i_node(:)  !! node-based work array

    ! Floating "pointers" for all arrays
    integer, private :: last_r_cell  !! pointer to last used real cell-array
    integer, private :: last_r_face  !! pointer to last used real face-array
    integer, private :: last_r_node  !! pointer to last used real node-array

    integer, private :: last_i_cell  !! pointer to last used integer cell-array
    integer, private :: last_i_face  !! pointer to last used integer face-array
    integer, private :: last_i_node  !! pointer to last used integer node-array

    ! Maximum number of cells, boundary cells, faces and nodes over all grids
    integer, private :: max_nc
    integer, private :: max_nb
    integer, private :: max_nf
    integer, private :: max_nn

    ! Maximum number of used arrays
    integer, private :: max_r_cell = 0  !! number of real cell-arrays used
    integer, private :: max_r_face = 0  !! number of real face-arrays used
    integer, private :: max_r_node = 0  !! number of real node-arrays used

    integer, private :: max_i_cell = 0  !! number of integer cell-arrays used
    integer, private :: max_i_face = 0  !! number of integer face-arrays used
    integer, private :: max_i_node = 0  !! number of integer node-arrays used

    contains
      procedure, private :: Allocate_Int_Cell
      procedure, private :: Allocate_Int_Face
      procedure, private :: Allocate_Int_Node
      procedure, private :: Allocate_Real_Cell
      procedure, private :: Allocate_Real_Face
      procedure, private :: Allocate_Real_Node
      procedure          :: Allocate_Work
      procedure          :: Connect_Int_Cell
      procedure          :: Connect_Int_Face
      procedure          :: Connect_Int_Node
      procedure          :: Connect_Real_Cell
      procedure          :: Connect_Real_Face
      procedure          :: Connect_Real_Node
      procedure          :: Safe_Connect_Int_Cell
      procedure          :: Safe_Connect_Int_Face
      procedure          :: Safe_Connect_Int_Node
      procedure          :: Safe_Connect_Real_Cell
      procedure          :: Safe_Connect_Real_Face
      procedure          :: Safe_Connect_Real_Node
      procedure          :: Unsafe_Connect_Int_Cell
      procedure          :: Unsafe_Connect_Int_Face
      procedure          :: Unsafe_Connect_Int_Node
      procedure          :: Unsafe_Connect_Real_Cell
      procedure          :: Unsafe_Connect_Real_Face
      procedure          :: Unsafe_Connect_Real_Node
      procedure          :: Disconnect_Int_Cell
      procedure          :: Disconnect_Int_Face
      procedure          :: Disconnect_Int_Node
      procedure          :: Disconnect_Real_Cell
      procedure          :: Disconnect_Real_Face
      procedure          :: Disconnect_Real_Node
      procedure          :: Finalize_Work

  end type

  ! Singleton object Work
  type(Work_Type) :: Work  !! singleton object Work introduced for seamless
    !! synchronization of local data and easy access to public member function

  contains

#   include "Work_Mod/Allocate_Int_Cell.f90"
#   include "Work_Mod/Allocate_Int_Face.f90"
#   include "Work_Mod/Allocate_Int_Node.f90"
#   include "Work_Mod/Allocate_Real_Cell.f90"
#   include "Work_Mod/Allocate_Real_Face.f90"
#   include "Work_Mod/Allocate_Real_Node.f90"
#   include "Work_Mod/Allocate_Work.f90"
#   include "Work_Mod/Connect_Int_Cell.f90"
#   include "Work_Mod/Connect_Int_Face.f90"
#   include "Work_Mod/Connect_Int_Node.f90"
#   include "Work_Mod/Connect_Real_Cell.f90"
#   include "Work_Mod/Connect_Real_Face.f90"
#   include "Work_Mod/Connect_Real_Node.f90"
#   include "Work_Mod/Safe/Connect_Int_Cell.f90"
#   include "Work_Mod/Safe/Connect_Int_Face.f90"
#   include "Work_Mod/Safe/Connect_Int_Node.f90"
#   include "Work_Mod/Safe/Connect_Real_Cell.f90"
#   include "Work_Mod/Safe/Connect_Real_Face.f90"
#   include "Work_Mod/Safe/Connect_Real_Node.f90"
#   include "Work_Mod/Unsafe/Connect_Int_Cell.f90"
#   include "Work_Mod/Unsafe/Connect_Int_Face.f90"
#   include "Work_Mod/Unsafe/Connect_Int_Node.f90"
#   include "Work_Mod/Unsafe/Connect_Real_Cell.f90"
#   include "Work_Mod/Unsafe/Connect_Real_Face.f90"
#   include "Work_Mod/Unsafe/Connect_Real_Node.f90"
#   include "Work_Mod/Disconnect_Int_Cell.f90"
#   include "Work_Mod/Disconnect_Int_Face.f90"
#   include "Work_Mod/Disconnect_Int_Node.f90"
#   include "Work_Mod/Disconnect_Real_Cell.f90"
#   include "Work_Mod/Disconnect_Real_Face.f90"
#   include "Work_Mod/Disconnect_Real_Node.f90"
#   include "Work_Mod/Finalize_Work.f90"

  end module
