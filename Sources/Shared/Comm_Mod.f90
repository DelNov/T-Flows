#include "../Shared/Unused.h90"

!==============================================================================!
  module Comm_Mod
!------------------------------------------------------------------------------!
!>  The Comm_Mod module in T-Flows is a framework for managing Message Passing
!>  Interface (MPI) functionality in parallel computing environments.
!>  It provides an array of subroutines and types to facilitate effective
!>  communication, data exchange and syncronizations between processors.
!>  One of the important components are data members and procedure for handling
!>  parallel I/O, used extensivelly (and exclusivelly) to manage backup files.
!>  The module seamlessly integrates both parallel and sequential executions,
!>  ensuring adaptability across various computing environments.
!>  Additionally, it includes buffer and communicator types for efficient
!>  handling of distributed data and control over parallel processes.
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
  !> Encapsulates data which facilitates management
  !> of communication between different processors.
  type Buffer_Type
    integer              :: n_items    !! number if items
    integer, allocatable :: map(:)     !! map to local items
    integer, allocatable :: i_buff(:)  !! integer values stored in buffers
    logical, allocatable :: l_buff(:)  !! logical values stored in buffers
    real,    allocatable :: r_buff(:)  !! real values stored in buffers
    real,    allocatable :: o_buff(:)  !! old real values stored in buffers
  end type

  !---------------!
  !   Comm type   !
  !---------------!
  !> Encapsulates all data and procedures related to MPI functionality.
  !> It manages various aspects of parallel computing, including processor
  !> identification, buffer management, and communication between processors.
  !> One instance of Comm_Type is defined globally (called Global) for global
  !> collective functions, but one copy is also defined inside the Grid object
  !> for exchanging data between processors and saving backup files.
  type Comm_Type

    ! Number of buffer cells
    integer :: n_buff_cells  !! number of buffer cells in the sub-domain

    ! Processor i.d. defined for each cell
    integer, allocatable :: cell_proc(:)  !! processor rank at each cell

    ! Global cell and node numbers
    integer, allocatable :: cell_glo(:)  !! global cell numbers
    integer, allocatable :: node_glo(:)  !! global node numbers

    ! Aiding the backup saving to single file
    integer :: nc_sub   !! number of cells in subdomain
    integer :: nb_sub   !! number of bundary cells in subdomain
    integer :: nb_f     !! first boundary cell to save
    integer :: nb_l     !! last boundary cell to save
    integer :: nc_tot   !! total number of cells
    integer :: nb_tot   !! total number of bundary cells

    ! Single precision coud not be avoided here :-(
    integer(SP), allocatable :: cell_map(:)      !! cell map for parallel I/O
    integer(SP), allocatable :: bnd_cell_map(:)  !! boundary cell map for
                                                 !! parallel I/O
#   if T_FLOWS_MPI == 1
      type(Mpi_Datatype), private :: cell_map_type
      type(Mpi_Datatype), private :: bnd_cell_map_type
#   else
      integer, private :: cell_map_type
      integer, private :: bnd_cell_map_type
#   endif

    ! Number of processors per node and processor i.d.s for each node
    type(Buffer_Type), allocatable :: cells_send(:)  !! send buffers
    type(Buffer_Type), allocatable :: cells_recv(:)  !! receive buffers

    integer, private :: n_processors    !! number of processors
    integer, private :: this_processor  !! current processor

    contains

      ! File management
      procedure :: Close_File
      procedure :: Create_New_Types
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
      procedure :: Exchange_Int_Array
      procedure :: Exchange_Log_Array
      procedure :: Exchange_Real_Array
      procedure :: Sendrecv_Int_Arrays
      procedure :: Sendrecv_Log_Arrays
      procedure :: Sendrecv_Real_Arrays

      ! Global
      procedure :: End_Parallel
      procedure :: Lor_Log
      procedure :: Lor_Log_Array
      procedure :: Max_Int
      procedure :: Max_Real
      procedure :: Max_Real_Array
      procedure :: Max_Reals
      procedure :: Min_Int
      procedure :: Min_Real
      procedure :: Start_Parallel
      procedure :: Sum_Int
      procedure :: Sum_Int_Array
      procedure :: Sum_Ints
      procedure :: Sum_Real
      procedure :: Sum_Real_Array
      procedure :: Sum_Reals
      procedure :: Wait
  end type

  !------------------------------------------------------------------------!
  !   A big global communicator, introduced essentiall to give access to   !
  !   private variables n_processors and this_processor to other objects   !
  !------------------------------------------------------------------------!
  type(Comm_Type) :: Global

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

    ! Shared pure function
#   include "Comm_Mod/Shared/First_Proc.f90"
#   include "Comm_Mod/Shared/N_Procs.f90"
#   include "Comm_Mod/Shared/Parallel_Run.f90"
#   include "Comm_Mod/Shared/Sequential_Run.f90"
#   include "Comm_Mod/Shared/This_Proc.f90"

# if T_FLOWS_MPI == 1
#   include "Comm_Mod/Parallel/Global/Start_Parallel.f90"
#   include "Comm_Mod/Parallel/Global/Wait.f90"
#   include "Comm_Mod/Parallel/Global/End_Parallel.f90"

    ! Parallel I/O (used for backup file management)
#   include "Comm_Mod/Parallel/Input_Output/Close_File.f90"
#   include "Comm_Mod/Parallel/Input_Output/Create_New_Types.f90"
#   include "Comm_Mod/Parallel/Input_Output/Open_File_Read.f90"
#   include "Comm_Mod/Parallel/Input_Output/Open_File_Write.f90"
#   include "Comm_Mod/Parallel/Input_Output/Read_Int.f90"
#   include "Comm_Mod/Parallel/Input_Output/Read_Int_Array.f90"
#   include "Comm_Mod/Parallel/Input_Output/Read_Log.f90"
#   include "Comm_Mod/Parallel/Input_Output/Read_Log_Array.f90"
#   include "Comm_Mod/Parallel/Input_Output/Read_Bnd_Real.f90"
#   include "Comm_Mod/Parallel/Input_Output/Read_Cell_Real.f90"
#   include "Comm_Mod/Parallel/Input_Output/Read_Real.f90"
#   include "Comm_Mod/Parallel/Input_Output/Read_Real_Array.f90"
#   include "Comm_Mod/Parallel/Input_Output/Read_Text.f90"
#   include "Comm_Mod/Parallel/Input_Output/Write_Int.f90"
#   include "Comm_Mod/Parallel/Input_Output/Write_Int_Array.f90"
#   include "Comm_Mod/Parallel/Input_Output/Write_Log.f90"
#   include "Comm_Mod/Parallel/Input_Output/Write_Log_Array.f90"
#   include "Comm_Mod/Parallel/Input_Output/Write_Bnd_Real.f90"
#   include "Comm_Mod/Parallel/Input_Output/Write_Cell_Real.f90"
#   include "Comm_Mod/Parallel/Input_Output/Write_Real.f90"
#   include "Comm_Mod/Parallel/Input_Output/Write_Real_Array.f90"
#   include "Comm_Mod/Parallel/Input_Output/Write_Text.f90"

    ! Global communicatins are better of as non-members
#   include "Comm_Mod/Parallel/Global/Lor_Log.f90"
#   include "Comm_Mod/Parallel/Global/Lor_Log_Array.f90"
#   include "Comm_Mod/Parallel/Global/Max_Real.f90"
#   include "Comm_Mod/Parallel/Global/Max_Real_Array.f90"
#   include "Comm_Mod/Parallel/Global/Max_Reals.f90"
#   include "Comm_Mod/Parallel/Global/Min_Real.f90"
#   include "Comm_Mod/Parallel/Global/Max_Int.f90"
#   include "Comm_Mod/Parallel/Global/Min_Int.f90"
#   include "Comm_Mod/Parallel/Global/Sum_Int.f90"
#   include "Comm_Mod/Parallel/Global/Sum_Int_Array.f90"
#   include "Comm_Mod/Parallel/Global/Sum_Ints.f90"
#   include "Comm_Mod/Parallel/Global/Sum_Real.f90"
#   include "Comm_Mod/Parallel/Global/Sum_Real_Array.f90"
#   include "Comm_Mod/Parallel/Global/Sum_Reals.f90"

    ! Messaging
#   include "Comm_Mod/Parallel/Messaging/Exchange_Int_Array.f90"
#   include "Comm_Mod/Parallel/Messaging/Exchange_Log_Array.f90"
#   include "Comm_Mod/Parallel/Messaging/Exchange_Real_Array.f90"
#   include "Comm_Mod/Parallel/Messaging/Sendrecv_Int_Arrays.f90"
#   include "Comm_Mod/Parallel/Messaging/Sendrecv_Log_Arrays.f90"
#   include "Comm_Mod/Parallel/Messaging/Sendrecv_Real_Arrays.f90"
# else
#   include "Comm_Mod/Sequential/Global/Start_Parallel.f90"
#   include "Comm_Mod/Sequential/Global/Wait.f90"
#   include "Comm_Mod/Sequential/Global/End_Parallel.f90"

    ! Sequential I/O (used for backup file management)
#   include "Comm_Mod/Sequential/Input_Output/Close_File.f90"
#   include "Comm_Mod/Sequential/Input_Output/Create_New_Types.f90"
#   include "Comm_Mod/Sequential/Input_Output/Open_File_Read.f90"
#   include "Comm_Mod/Sequential/Input_Output/Open_File_Write.f90"
#   include "Comm_Mod/Sequential/Input_Output/Read_Int.f90"
#   include "Comm_Mod/Sequential/Input_Output/Read_Int_Array.f90"
#   include "Comm_Mod/Sequential/Input_Output/Read_Log.f90"
#   include "Comm_Mod/Sequential/Input_Output/Read_Log_Array.f90"
#   include "Comm_Mod/Sequential/Input_Output/Read_Bnd_Real.f90"
#   include "Comm_Mod/Sequential/Input_Output/Read_Cell_Real.f90"
#   include "Comm_Mod/Sequential/Input_Output/Read_Real.f90"
#   include "Comm_Mod/Sequential/Input_Output/Read_Real_Array.f90"
#   include "Comm_Mod/Sequential/Input_Output/Read_Text.f90"
#   include "Comm_Mod/Sequential/Input_Output/Write_Int.f90"
#   include "Comm_Mod/Sequential/Input_Output/Write_Int_Array.f90"
#   include "Comm_Mod/Sequential/Input_Output/Write_Log.f90"
#   include "Comm_Mod/Sequential/Input_Output/Write_Log_Array.f90"
#   include "Comm_Mod/Sequential/Input_Output/Write_Bnd_Real.f90"
#   include "Comm_Mod/Sequential/Input_Output/Write_Cell_Real.f90"
#   include "Comm_Mod/Sequential/Input_Output/Write_Real.f90"
#   include "Comm_Mod/Sequential/Input_Output/Write_Real_Array.f90"
#   include "Comm_Mod/Sequential/Input_Output/Write_Text.f90"

    ! Global communicatins are better of as non-members
#   include "Comm_Mod/Sequential/Global/Lor_Log.f90"
#   include "Comm_Mod/Sequential/Global/Lor_Log_Array.f90"
#   include "Comm_Mod/Sequential/Global/Max_Real.f90"
#   include "Comm_Mod/Sequential/Global/Max_Real_Array.f90"
#   include "Comm_Mod/Sequential/Global/Max_Reals.f90"
#   include "Comm_Mod/Sequential/Global/Min_Real.f90"
#   include "Comm_Mod/Sequential/Global/Max_Int.f90"
#   include "Comm_Mod/Sequential/Global/Min_Int.f90"
#   include "Comm_Mod/Sequential/Global/Sum_Int.f90"
#   include "Comm_Mod/Sequential/Global/Sum_Int_Array.f90"
#   include "Comm_Mod/Sequential/Global/Sum_Ints.f90"
#   include "Comm_Mod/Sequential/Global/Sum_Real.f90"
#   include "Comm_Mod/Sequential/Global/Sum_Real_Array.f90"
#   include "Comm_Mod/Sequential/Global/Sum_Reals.f90"

    ! Messaging
#   include "Comm_Mod/Sequential/Messaging/Exchange_Int_Array.f90"
#   include "Comm_Mod/Sequential/Messaging/Exchange_Log_Array.f90"
#   include "Comm_Mod/Sequential/Messaging/Exchange_Real_Array.f90"
#   include "Comm_Mod/Sequential/Messaging/Sendrecv_Int_Arrays.f90"
#   include "Comm_Mod/Sequential/Messaging/Sendrecv_Log_Arrays.f90"
#   include "Comm_Mod/Sequential/Messaging/Sendrecv_Real_Arrays.f90"
# endif

  end module
