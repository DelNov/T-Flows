#include "../Shared/Assert.h90"

#define CHECK_USAGE 1

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

  !---------------------!
  !  Real poiner type   !
  !---------------------!
  !> Simple container of a real, contiguous pointer.  It was introduced to
  !> facilitate creation of variable size arrays containing pointers only.
  type Real_Pointer_Type
    real, pointer, contiguous :: ptr(:)
  end type

  !------------------------!
  !  Integer poiner type   !
  !------------------------!
  !> Simple container of an integer contiguous pointer.  It was introduced to
  !> facilitate creation of variable size arrays containing pointers only.
  type Int_Pointer_Type
    integer, pointer, contiguous :: ptr(:)
  end type

  !--------------!
  !  Work type   !
  !--------------!
  !> The Work_Type within this module encapsulates various arrays for both
  !> integer and real data types, specifically tailored for different grid
  !> elements like cells, faces, and nodes.
  type Work_Type

    ! Real working arrays
    type(Real_Pointer_Type), allocatable :: r_cell(:)  !! cell-based work array
    type(Real_Pointer_Type), allocatable :: r_face(:)  !! face-based work array
    type(Real_Pointer_Type), allocatable :: r_node(:)  !! node-based work array

    ! Integer working arrays
    type(Int_Pointer_Type), allocatable :: i_cell(:)  !! cell-based work array
    type(Int_Pointer_Type), allocatable :: i_face(:)  !! face-based work array
    type(Int_Pointer_Type), allocatable :: i_node(:)  !! node-based work array

    ! Floating "pointers" for all arrays
    integer, private :: last_r_cell  !! pointer to last used real cell-array
    integer, private :: last_r_face  !! pointer to last used real face-array
    integer, private :: last_r_node  !! pointer to last used real node-array

    integer, private :: last_i_cell  !! pointer to last used integer cell-array
    integer, private :: last_i_face  !! pointer to last used integer face-array
    integer, private :: last_i_node  !! pointer to last used integer node-array

    ! Maximum number of used arrays
#ifdef CHECK_USAGE
    integer, private :: max_r_cell = 0  !! number of real cell-arrays used
    integer, private :: max_r_face = 0  !! number of real face-arrays used
    integer, private :: max_r_node = 0  !! number of real node-arrays used

    integer, private :: max_i_cell = 0  !! number of integer cell-arrays used
    integer, private :: max_i_face = 0  !! number of integer face-arrays used
    integer, private :: max_i_node = 0  !! number of integer node-arrays used
#endif

    ! Requested number of used arrays
    integer, private :: req_r_cell = 0  !! number of requested real cell-arrays
    integer, private :: req_r_face = 0  !! number of requested real face-arrays
    integer, private :: req_r_node = 0  !! number of requested real node-arrays

    integer, private :: req_i_cell = 0  !! number of requested int. cell-arrays
    integer, private :: req_i_face = 0  !! number of requested int. face-arrays
    integer, private :: req_i_node = 0  !! number of requested int. node-arrays

    logical :: allocated = .false.

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
#ifdef CHECK_USAGE
#   include "Work_Mod/Check_Usage/Connect_Int_Cell.f90"
#   include "Work_Mod/Check_Usage/Connect_Int_Face.f90"
#   include "Work_Mod/Check_Usage/Connect_Int_Node.f90"
#   include "Work_Mod/Check_Usage/Connect_Real_Cell.f90"
#   include "Work_Mod/Check_Usage/Connect_Real_Face.f90"
#   include "Work_Mod/Check_Usage/Connect_Real_Node.f90"
#   include "Work_Mod/Check_Usage/Disconnect_Int_Cell.f90"
#   include "Work_Mod/Check_Usage/Disconnect_Int_Face.f90"
#   include "Work_Mod/Check_Usage/Disconnect_Int_Node.f90"
#   include "Work_Mod/Check_Usage/Disconnect_Real_Cell.f90"
#   include "Work_Mod/Check_Usage/Disconnect_Real_Face.f90"
#   include "Work_Mod/Check_Usage/Disconnect_Real_Node.f90"
#   include "Work_Mod/Check_Usage/Finalize_Work.f90"
#else
#   include "Work_Mod/No_Checking/Connect_Int_Cell.f90"
#   include "Work_Mod/No_Checking/Connect_Int_Face.f90"
#   include "Work_Mod/No_Checking/Connect_Int_Node.f90"
#   include "Work_Mod/No_Checking/Connect_Real_Cell.f90"
#   include "Work_Mod/No_Checking/Connect_Real_Face.f90"
#   include "Work_Mod/No_Checking/Connect_Real_Node.f90"
#   include "Work_Mod/No_Checking/Disconnect_Int_Cell.f90"
#   include "Work_Mod/No_Checking/Disconnect_Int_Face.f90"
#   include "Work_Mod/No_Checking/Disconnect_Int_Node.f90"
#   include "Work_Mod/No_Checking/Disconnect_Real_Cell.f90"
#   include "Work_Mod/No_Checking/Disconnect_Real_Face.f90"
#   include "Work_Mod/No_Checking/Disconnect_Real_Node.f90"
#   include "Work_Mod/No_Checking/Finalize_Work.f90"
#endif

  end module
