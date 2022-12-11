!==============================================================================!
  module Comm_Mod
!------------------------------------------------------------------------------!
!   Module for MPI functionality.                                              !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
# if T_FLOWS_MPI == 1
    use Mpi_f08
# endif
  use Const_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !-----------------!
  !   Buffer type   !
  !-----------------!
  type Buffer_Type
    integer              :: n_items
    integer, allocatable :: map(:)     ! map to local items
    integer, allocatable :: i_buff(:)  ! integer values stored in buffers
    logical, allocatable :: l_buff(:)  ! logical values stored in buffers
    real,    allocatable :: r_buff(:)  ! real values stored in buffers
    real,    allocatable :: o_buff(:)  ! old real values stored in buffers
  end type

  !---------------!
  !   Comm type   !
  !---------------!
  type Comm_Type    ! also used inside the Grid_Type

    ! Number of buffer cells
    integer :: n_buff_cells

    ! Processor i.d. defined for each cell
    integer, allocatable :: cell_proc(:)

    ! Global cell and node numbers
    integer, contiguous, pointer :: cell_glo(:)
    integer, contiguous, pointer :: node_glo(:)

    ! Variables which follow are for backup saving to single file
    integer :: nc_sub   ! number of cells in subdomain
    integer :: nb_sub   ! number of bundary cells in subdomain
    integer :: nb_f     ! first boundary cell to save
    integer :: nb_l     ! last boundary cell to save
    integer :: nc_tot   ! total number of cells
    integer :: nb_tot   ! total number of bundary cells

    ! Single precision coud not be avoided here :-(
    integer(SP), allocatable :: cell_map(:)
    integer(SP), allocatable :: bnd_cell_map(:)

#   if T_FLOWS_MPI == 1
      type(Mpi_Datatype), private :: cell_map_type
      type(Mpi_Datatype), private :: bnd_cell_map_type
#   else
      integer, private :: cell_map_type
      integer, private :: bnd_cell_map_type
#   endif

    ! Number of processors per node and processor i.d.s for each node
    type(Buffer_Type), allocatable :: cells_send(:)
    type(Buffer_Type), allocatable :: cells_recv(:)

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

  integer :: this_proc  ! processor i.d.
  integer :: n_proc     ! number of processors

  ! These communication types will depend on precision
#if T_FLOWS_MPI == 1
  type(Mpi_Datatype) :: comm_type_int
  type(Mpi_Datatype) :: comm_type_log
  type(Mpi_Datatype) :: comm_type_real
#else
  integer :: comm_type_int
  integer :: comm_type_log
  integer :: comm_type_real
#endif

  contains

# if T_FLOWS_MPI == 1
    ! Three basic ones are non-member
#   include "Comm_Mod/Parallel/Start.f90"
#   include "Comm_Mod/Parallel/Wait.f90"
#   include "Comm_Mod/Parallel/End.f90"

    ! File management
#   include "Comm_Mod/Parallel/Close_File.f90"
#   include "Comm_Mod/Parallel/Open_File_Read.f90"
#   include "Comm_Mod/Parallel/Open_File_Write.f90"
#   include "Comm_Mod/Parallel/Read_Int.f90"
#   include "Comm_Mod/Parallel/Read_Int_Array.f90"
#   include "Comm_Mod/Parallel/Read_Log.f90"
#   include "Comm_Mod/Parallel/Read_Log_Array.f90"
#   include "Comm_Mod/Parallel/Read_Bnd_Real.f90"
#   include "Comm_Mod/Parallel/Read_Cell_Real.f90"
#   include "Comm_Mod/Parallel/Read_Real.f90"
#   include "Comm_Mod/Parallel/Read_Real_Array.f90"
#   include "Comm_Mod/Parallel/Read_Text.f90"
#   include "Comm_Mod/Parallel/Write_Int.f90"
#   include "Comm_Mod/Parallel/Write_Int_Array.f90"
#   include "Comm_Mod/Parallel/Write_Log.f90"
#   include "Comm_Mod/Parallel/Write_Log_Array.f90"
#   include "Comm_Mod/Parallel/Write_Bnd_Real.f90"
#   include "Comm_Mod/Parallel/Write_Cell_Real.f90"
#   include "Comm_Mod/Parallel/Write_Real.f90"
#   include "Comm_Mod/Parallel/Write_Real_Array.f90"
#   include "Comm_Mod/Parallel/Write_Text.f90"

    ! Global communicatins are better of as non-members
#   include "Comm_Mod/Parallel/Global_Lor_Log.f90"
#   include "Comm_Mod/Parallel/Global_Lor_Log_Array.f90"
#   include "Comm_Mod/Parallel/Global_Max_Real.f90"
#   include "Comm_Mod/Parallel/Global_Min_Real.f90"
#   include "Comm_Mod/Parallel/Global_Max_Int.f90"
#   include "Comm_Mod/Parallel/Global_Min_Int.f90"
#   include "Comm_Mod/Parallel/Global_Sum_Int.f90"
#   include "Comm_Mod/Parallel/Global_Sum_Int_Array.f90"
#   include "Comm_Mod/Parallel/Global_Sum_Real.f90"
#   include "Comm_Mod/Parallel/Global_Sum_Real_Array.f90"

    ! Messaging
#   include "Comm_Mod/Parallel/Create_New_Types.f90"
#   include "Comm_Mod/Parallel/Exchange_Int_Array.f90"
#   include "Comm_Mod/Parallel/Exchange_Log_Array.f90"
#   include "Comm_Mod/Parallel/Exchange_Real_Array.f90"
#   include "Comm_Mod/Parallel/Sendrecv_Int_Arrays.f90"
#   include "Comm_Mod/Parallel/Sendrecv_Log_Arrays.f90"
#   include "Comm_Mod/Parallel/Sendrecv_Real_Arrays.f90"
# else
    ! Three basic ones are non-member
#   include "Comm_Mod/Sequential/Start.f90"
#   include "Comm_Mod/Sequential/Wait.f90"
#   include "Comm_Mod/Sequential/End.f90"

    ! File management
#   include "Comm_Mod/Sequential/Close_File.f90"
#   include "Comm_Mod/Sequential/Open_File_Read.f90"
#   include "Comm_Mod/Sequential/Open_File_Write.f90"
#   include "Comm_Mod/Sequential/Read_Int.f90"
#   include "Comm_Mod/Sequential/Read_Int_Array.f90"
#   include "Comm_Mod/Sequential/Read_Log.f90"
#   include "Comm_Mod/Sequential/Read_Log_Array.f90"
#   include "Comm_Mod/Sequential/Read_Bnd_Real.f90"
#   include "Comm_Mod/Sequential/Read_Cell_Real.f90"
#   include "Comm_Mod/Sequential/Read_Real.f90"
#   include "Comm_Mod/Sequential/Read_Real_Array.f90"
#   include "Comm_Mod/Sequential/Read_Text.f90"
#   include "Comm_Mod/Sequential/Write_Int.f90"
#   include "Comm_Mod/Sequential/Write_Int_Array.f90"
#   include "Comm_Mod/Sequential/Write_Log.f90"
#   include "Comm_Mod/Sequential/Write_Log_Array.f90"
#   include "Comm_Mod/Sequential/Write_Bnd_Real.f90"
#   include "Comm_Mod/Sequential/Write_Cell_Real.f90"
#   include "Comm_Mod/Sequential/Write_Real.f90"
#   include "Comm_Mod/Sequential/Write_Real_Array.f90"
#   include "Comm_Mod/Sequential/Write_Text.f90"

    ! Global communicatins are better of as non-members
#   include "Comm_Mod/Sequential/Global_Lor_Log.f90"
#   include "Comm_Mod/Sequential/Global_Lor_Log_Array.f90"
#   include "Comm_Mod/Sequential/Global_Max_Real.f90"
#   include "Comm_Mod/Sequential/Global_Min_Real.f90"
#   include "Comm_Mod/Sequential/Global_Max_Int.f90"
#   include "Comm_Mod/Sequential/Global_Min_Int.f90"
#   include "Comm_Mod/Sequential/Global_Sum_Int.f90"
#   include "Comm_Mod/Sequential/Global_Sum_Int_Array.f90"
#   include "Comm_Mod/Sequential/Global_Sum_Real.f90"
#   include "Comm_Mod/Sequential/Global_Sum_Real_Array.f90"

    ! Messaging
#   include "Comm_Mod/Sequential/Create_New_Types.f90"
#   include "Comm_Mod/Sequential/Exchange_Int_Array.f90"
#   include "Comm_Mod/Sequential/Exchange_Log_Array.f90"
#   include "Comm_Mod/Sequential/Exchange_Real_Array.f90"
#   include "Comm_Mod/Sequential/Sendrecv_Int_Arrays.f90"
#   include "Comm_Mod/Sequential/Sendrecv_Log_Arrays.f90"
#   include "Comm_Mod/Sequential/Sendrecv_Real_Arrays.f90"
# endif

  end module
