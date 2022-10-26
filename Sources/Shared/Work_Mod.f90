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

  !---------------------!
  !  Real poiner type   !
  !---------------------!
  type Real_Pointer_Type
    real, pointer, contiguous :: ptr(:)
  end type

  !------------------------!
  !  Integer poiner type   !
  !------------------------!
  type Int_Pointer_Type
    integer, pointer, contiguous :: ptr(:)
  end type

  !--------------!
  !  Work type   !
  !--------------!
  type Work_Type

    ! Real working arrays
    type(Real_Pointer_Type), allocatable :: r_cell(:)
    type(Real_Pointer_Type), allocatable :: r_face(:)
    type(Real_Pointer_Type), allocatable :: r_node(:)

    ! Integer working arrays
    type(Int_Pointer_Type), allocatable :: i_cell(:)
    type(Int_Pointer_Type), allocatable :: i_face(:)
    type(Int_Pointer_Type), allocatable :: i_node(:)

    integer :: last_r_cell
    integer :: last_r_face
    integer :: last_r_node

    integer :: last_i_cell
    integer :: last_i_face
    integer :: last_i_node

    integer :: max_r_cell = 0
    integer :: max_r_face = 0
    integer :: max_r_node = 0

    integer :: max_i_cell = 0
    integer :: max_i_face = 0
    integer :: max_i_node = 0

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

  end type

  type(Work_Type) :: Work

  contains

  include 'Work_Mod/Allocate_Int_Cell.f90'
  include 'Work_Mod/Allocate_Int_Face.f90'
  include 'Work_Mod/Allocate_Int_Node.f90'
  include 'Work_Mod/Allocate_Real_Cell.f90'
  include 'Work_Mod/Allocate_Real_Face.f90'
  include 'Work_Mod/Allocate_Real_Node.f90'
  include 'Work_Mod/Allocate_Work.f90'
  include 'Work_Mod/Connect_Int_Cell.f90'
  include 'Work_Mod/Connect_Int_Face.f90'
  include 'Work_Mod/Connect_Int_Node.f90'
  include 'Work_Mod/Connect_Real_Cell.f90'
  include 'Work_Mod/Connect_Real_Face.f90'
  include 'Work_Mod/Connect_Real_Node.f90'
  include 'Work_Mod/Disconnect_Int_Cell.f90'
  include 'Work_Mod/Disconnect_Int_Face.f90'
  include 'Work_Mod/Disconnect_Int_Node.f90'
  include 'Work_Mod/Disconnect_Real_Cell.f90'
  include 'Work_Mod/Disconnect_Real_Face.f90'
  include 'Work_Mod/Disconnect_Real_Node.f90'

  end module
