!==============================================================================!
  module Comm_Mod
!------------------------------------------------------------------------------!
!   Module for MPI functionality.                                              !
!------------------------------------------------------------------------------!
  implicit none
!----------------------------------[Include]-----------------------------------!
  include 'mpif.h'
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
    real,    allocatable :: r_buff(:)  ! real    values stored in buffers
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

    ! (kind=4) coud not be avoided here :-(
    integer(kind=4), allocatable :: cell_map(:)
    integer(kind=4), allocatable :: bnd_cell_map(:)

    ! Variables which follow are for backup saving to single file
    integer :: nc_s   ! number of cells in subdomain
    integer :: nb_s   ! number of bundary cells in subdomain
    integer :: nb_f   ! first boundary cell to save
    integer :: nb_l   ! last boundary cell to save
    integer :: nc_t   ! total number of cells 
    integer :: nb_t   ! total number of bundary cells

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

  integer :: cell_map_type
  integer :: bnd_cell_map_type

  contains

  include 'Comm_Mod/Parallel/Close_File.f90'
  include 'Comm_Mod/Parallel/Create_New_Types.f90'
  include 'Comm_Mod/Parallel/End.f90'
  include 'Comm_Mod/Parallel/Exchange_Int_Array.f90'
  include 'Comm_Mod/Parallel/Exchange_Log_Array.f90'
  include 'Comm_Mod/Parallel/Exchange_Real_Array.f90'
  include 'Comm_Mod/Parallel/Global_Lor_Log_Array.f90'
  include 'Comm_Mod/Parallel/Global_Max_Real.f90'
  include 'Comm_Mod/Parallel/Global_Min_Real.f90'
  include 'Comm_Mod/Parallel/Global_Max_Int.f90'
  include 'Comm_Mod/Parallel/Global_Min_Int.f90'
  include 'Comm_Mod/Parallel/Global_Sum_Int.f90'
  include 'Comm_Mod/Parallel/Global_Sum_Int_Array.f90'
  include 'Comm_Mod/Parallel/Global_Sum_Real.f90'
  include 'Comm_Mod/Parallel/Global_Sum_Real_Array.f90'
  include 'Comm_Mod/Parallel/Open_File_Read.f90'
  include 'Comm_Mod/Parallel/Open_File_Write.f90'
  include 'Comm_Mod/Parallel/Read_Int.f90'
  include 'Comm_Mod/Parallel/Read_Int_Array.f90'
  include 'Comm_Mod/Parallel/Read_Log_Array.f90'
  include 'Comm_Mod/Parallel/Read_Bnd_Real.f90'
  include 'Comm_Mod/Parallel/Read_Cell_Real.f90'
  include 'Comm_Mod/Parallel/Read_Real.f90'
  include 'Comm_Mod/Parallel/Read_Real_Array.f90'
  include 'Comm_Mod/Parallel/Read_Text.f90'
  include 'Comm_Mod/Parallel/Recv_Int_Array.f90'
  include 'Comm_Mod/Parallel/Recv_Log_Array.f90'
  include 'Comm_Mod/Parallel/Recv_Real_Array.f90'
  include 'Comm_Mod/Parallel/Send_Int_Array.f90'
  include 'Comm_Mod/Parallel/Send_Log_Array.f90'
  include 'Comm_Mod/Parallel/Send_Real_Array.f90'
  include 'Comm_Mod/Parallel/Sendrecv_Int_Arrays.f90'
  include 'Comm_Mod/Parallel/Sendrecv_Log_Arrays.f90'
  include 'Comm_Mod/Parallel/Sendrecv_Real_Arrays.f90'
  include 'Comm_Mod/Parallel/Start.f90'
  include 'Comm_Mod/Parallel/Wait.f90'
  include 'Comm_Mod/Parallel/Write_Int.f90'
  include 'Comm_Mod/Parallel/Write_Int_Array.f90'
  include 'Comm_Mod/Parallel/Write_Log_Array.f90'
  include 'Comm_Mod/Parallel/Write_Bnd_Real.f90'
  include 'Comm_Mod/Parallel/Write_Cell_Real.f90'
  include 'Comm_Mod/Parallel/Write_Real.f90'
  include 'Comm_Mod/Parallel/Write_Real_Array.f90'
  include 'Comm_Mod/Parallel/Write_Text.f90'

  end module
