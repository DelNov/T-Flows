!==============================================================================!
  module Comm_Mod
!------------------------------------------------------------------------------!
!   Module for no MPI functionality.                                           !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
!------------------------------------------------------------------------------!
  implicit none
!----------------------------------[Include]-----------------------------------!
! include 'mpif.h'
!==============================================================================!

  !-----------------!
  !                 !
  !   Buffer type   !
  !                 !
  !-----------------!
  type Buffer_Type
    integer              :: n_items
    integer, allocatable :: map  (:)   ! map to local node / face
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
    integer(SP), allocatable :: cell_map(:)
    integer(SP), allocatable :: bnd_cell_map(:)
    integer(SP), allocatable :: face_map(:)

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

  end type

  ! Parameters for size of typical variables in bytes
  integer, parameter :: SIZE_INT  = 8
  integer, parameter :: SIZE_LOG  = 8
  integer, parameter :: SIZE_REAL = 8

  integer :: this_proc  ! processor i.d.
  integer :: n_proc     ! number of processors

  contains

  include 'Comm_Mod/Sequential/Close_File.f90'
  include 'Comm_Mod/Sequential/Create_New_Types.f90'
  include 'Comm_Mod/Sequential/End.f90'
  include 'Comm_Mod/Sequential/Exchange_Int_Array.f90'
  include 'Comm_Mod/Sequential/Exchange_Log_Array.f90'
  include 'Comm_Mod/Sequential/Exchange_Real_Array.f90'
  include 'Comm_Mod/Sequential/Global_Lor_Log_Array.f90'
  include 'Comm_Mod/Sequential/Global_Max_Real.f90'
  include 'Comm_Mod/Sequential/Global_Min_Real.f90'
  include 'Comm_Mod/Sequential/Global_Max_Int.f90'
  include 'Comm_Mod/Sequential/Global_Min_Int.f90'
  include 'Comm_Mod/Sequential/Global_Sum_Int.f90'
  include 'Comm_Mod/Sequential/Global_Sum_Int_Array.f90'
  include 'Comm_Mod/Sequential/Global_Sum_Real.f90'
  include 'Comm_Mod/Sequential/Global_Sum_Real_Array.f90'
  include 'Comm_Mod/Sequential/Open_File_Read.f90'
  include 'Comm_Mod/Sequential/Open_File_Write.f90'
  include 'Comm_Mod/Sequential/Read_Int.f90'
  include 'Comm_Mod/Sequential/Read_Int_Array.f90'
  include 'Comm_Mod/Sequential/Read_Log_Array.f90'
  include 'Comm_Mod/Sequential/Read_Bnd_Real.f90'
  include 'Comm_Mod/Sequential/Read_Cell_Real.f90'
  include 'Comm_Mod/Sequential/Read_Face_Real.f90'
  include 'Comm_Mod/Sequential/Read_Real.f90'
  include 'Comm_Mod/Sequential/Read_Real_Array.f90'
  include 'Comm_Mod/Sequential/Read_Text.f90'
  include 'Comm_Mod/Sequential/Recv_Int_Array.f90'
  include 'Comm_Mod/Sequential/Recv_Log_Array.f90'
  include 'Comm_Mod/Sequential/Recv_Real_Array.f90'
  include 'Comm_Mod/Sequential/Send_Int_Array.f90'
  include 'Comm_Mod/Sequential/Send_Log_Array.f90'
  include 'Comm_Mod/Sequential/Send_Real_Array.f90'
  include 'Comm_Mod/Sequential/Sendrecv_Int_Arrays.f90'
  include 'Comm_Mod/Sequential/Sendrecv_Log_Arrays.f90'
  include 'Comm_Mod/Sequential/Sendrecv_Real_Arrays.f90'
  include 'Comm_Mod/Sequential/Start.f90'
  include 'Comm_Mod/Sequential/Wait.f90'
  include 'Comm_Mod/Sequential/Write_Int.f90'
  include 'Comm_Mod/Sequential/Write_Int_Array.f90'
  include 'Comm_Mod/Sequential/Write_Log_Array.f90'
  include 'Comm_Mod/Sequential/Write_Bnd_Real.f90'
  include 'Comm_Mod/Sequential/Write_Cell_Real.f90'
  include 'Comm_Mod/Sequential/Write_Face_Real.f90'
  include 'Comm_Mod/Sequential/Write_Real.f90'
  include 'Comm_Mod/Sequential/Write_Real_Array.f90'
  include 'Comm_Mod/Sequential/Write_Text.f90'

  end module
