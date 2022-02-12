!==============================================================================!
  module Comm_Mod
!------------------------------------------------------------------------------!
!   Module for no MPI functionality.                                           !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
! use Mpi
  use Const_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !-----------------!
  !                 !
  !   Buffer type   !
  !                 !
  !-----------------!
  type Buffer_Type
    integer              :: n_items
    integer, allocatable :: map(:)     ! map to local node / face
    integer, allocatable :: i_buff(:)  ! integer values stored in buffers
    logical, allocatable :: l_buff(:)  ! logical values stored in buffers
    real,    allocatable :: r_buff(:)  ! real values stored in buffers
    real,    allocatable :: o_buff(:)  ! old real values stored in buffers
  end type

  !---------------!
  !               !
  !   Comm type   !
  !               !
  !---------------!
  type Comm_Type    ! also used inside the Grid_Type)

    ! Number of buffer cells
    integer :: n_buff_cells

    ! Processor i.d. defined for each cell
    integer, allocatable :: cell_proc(:)

    ! Global cell and node numbers
    integer, allocatable :: cell_glo(:)
    integer, allocatable :: node_glo(:)
    integer, allocatable :: face_glo(:)

    ! Single precision coud not be avoided here :-(
    integer, allocatable :: cell_map(:)
    integer, allocatable :: bnd_cell_map(:)

    ! Face maps.  They are not used in MPI calls, can be integers
    integer, allocatable :: face_map_uni_glo(:)  ! unique entries only
    integer, allocatable :: face_map_uni_loc(:)  ! unique entries only
    integer, allocatable :: face_map_dup_glo(:)  ! duplicate entries too
    integer, allocatable :: face_map_dup_loc(:)  ! duplicate entries too

    ! Variables which follow are for backup saving to single file
    integer :: nc_sub   ! number of cells in subdomain
    integer :: nb_sub   ! number of bundary cells in subdomain
    integer :: nb_f     ! first boundary cell to save
    integer :: nb_l     ! last boundary cell to save
    integer :: nc_tot   ! total number of cells
    integer :: nb_tot   ! total number of bundary cells
    integer :: nf_sub   ! number of faces in subdomain
    integer :: nf_tot   ! total number of faces

    ! Number of processors per node and processor i.d.s for each node
    integer,           allocatable :: nodes_n_procs(:)
    integer,           allocatable :: nodes_p(:,:)
    type(Buffer_Type), allocatable :: nodes_repl(:)
    type(Buffer_Type), allocatable :: cells_send(:)
    type(Buffer_Type), allocatable :: cells_recv(:)

    integer, private :: cell_map_type
    integer, private :: bnd_cell_map_type

    contains

      ! File management
      procedure :: Close_File
      procedure :: Open_File_Read
      procedure :: Open_File_Write
      procedure :: Read_Int
      procedure :: Read_Int_Array
      procedure :: Read_Log
      procedure :: Read_Log_Array
      procedure :: Read_Bnd_Real
      procedure :: Read_Cell_Real
      procedure :: Read_Real
      procedure :: Read_Real_Array
      procedure :: Read_Text
      procedure :: Write_Int
      procedure :: Write_Int_Array
      procedure :: Write_Log
      procedure :: Write_Log_Array
      procedure :: Write_Bnd_Real
      procedure :: Write_Cell_Real
      procedure :: Write_Real
      procedure :: Write_Real_Array
      procedure :: Write_Text

      ! Messaging
      procedure :: Create_New_Types
      procedure :: Exchange_Int_Array
      procedure :: Exchange_Log_Array
      procedure :: Exchange_Real_Array
      procedure :: Sendrecv_Int_Arrays
      procedure :: Sendrecv_Log_Arrays
      procedure :: Sendrecv_Real_Arrays

  end type

  ! Parameters for size of typical variables in bytes
  integer, parameter :: SIZE_INT  = 4
  integer, parameter :: SIZE_LOG  = 4
  integer, parameter :: SIZE_REAL = 8

  integer :: this_proc  ! processor i.d.
  integer :: n_proc     ! number of processors

  contains

  ! Three basic ones are non-member
  include 'Comm_Mod/Sequential/Start.f90'
  include 'Comm_Mod/Sequential/Wait.f90'
  include 'Comm_Mod/Sequential/End.f90'

  ! File management
  include 'Comm_Mod/Sequential/Close_File.f90'
  include 'Comm_Mod/Sequential/Open_File_Read.f90'
  include 'Comm_Mod/Sequential/Open_File_Write.f90'
  include 'Comm_Mod/Sequential/Read_Int.f90'
  include 'Comm_Mod/Sequential/Read_Int_Array.f90'
  include 'Comm_Mod/Sequential/Read_Log.f90'
  include 'Comm_Mod/Sequential/Read_Log_Array.f90'
  include 'Comm_Mod/Sequential/Read_Bnd_Real.f90'
  include 'Comm_Mod/Sequential/Read_Cell_Real.f90'
  include 'Comm_Mod/Sequential/Read_Real.f90'
  include 'Comm_Mod/Sequential/Read_Real_Array.f90'
  include 'Comm_Mod/Sequential/Read_Text.f90'
  include 'Comm_Mod/Sequential/Write_Int.f90'
  include 'Comm_Mod/Sequential/Write_Int_Array.f90'
  include 'Comm_Mod/Sequential/Write_Log.f90'
  include 'Comm_Mod/Sequential/Write_Log_Array.f90'
  include 'Comm_Mod/Sequential/Write_Bnd_Real.f90'
  include 'Comm_Mod/Sequential/Write_Cell_Real.f90'
  include 'Comm_Mod/Sequential/Write_Real.f90'
  include 'Comm_Mod/Sequential/Write_Real_Array.f90'
  include 'Comm_Mod/Sequential/Write_Text.f90'

  ! Global communicatins are better of as non-members
  include 'Comm_Mod/Sequential/Global_Lor_Log_Array.f90'
  include 'Comm_Mod/Sequential/Global_Max_Real.f90'
  include 'Comm_Mod/Sequential/Global_Min_Real.f90'
  include 'Comm_Mod/Sequential/Global_Max_Int.f90'
  include 'Comm_Mod/Sequential/Global_Min_Int.f90'
  include 'Comm_Mod/Sequential/Global_Sum_Int.f90'
  include 'Comm_Mod/Sequential/Global_Sum_Int_Array.f90'
  include 'Comm_Mod/Sequential/Global_Sum_Real.f90'
  include 'Comm_Mod/Sequential/Global_Sum_Real_Array.f90'

  ! Messaging
  include 'Comm_Mod/Sequential/Create_New_Types.f90'
  include 'Comm_Mod/Sequential/Exchange_Int_Array.f90'
  include 'Comm_Mod/Sequential/Exchange_Log_Array.f90'
  include 'Comm_Mod/Sequential/Exchange_Real_Array.f90'
  include 'Comm_Mod/Sequential/Sendrecv_Int_Arrays.f90'
  include 'Comm_Mod/Sequential/Sendrecv_Log_Arrays.f90'
  include 'Comm_Mod/Sequential/Sendrecv_Real_Arrays.f90'

  end module
