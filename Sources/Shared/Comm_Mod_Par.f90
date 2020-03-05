!==============================================================================!
  module Comm_Mod
!------------------------------------------------------------------------------!
!   Module for MPI functionality.                                              !
!------------------------------------------------------------------------------!
  implicit none
!----------------------------------[Include]-----------------------------------!
  include 'mpif.h'
!==============================================================================!

  !---------------!
  !               !
  !   Comm type   !
  !               !
  !---------------!
  type Comm_Type    ! also used inside the Grid_Type) 

    ! Number of buffer cells
    integer :: n_buff_cells

    ! For each buffer region: index of start (_s) and end (_e) buffer cell
    ! (Follows nomenclature in "../Shared/Bnd_Cond_Mod.f90")
    integer, allocatable :: buff_s_cell(:), buff_e_cell(:)

    ! Processor i.d. defined for each cell
    integer, allocatable :: cell_proc(:)

    ! Buffer index
    integer, allocatable :: buff_index(:)

    ! Global cell numbers
    integer, allocatable :: cell_glo(:)

    ! (kind=4) coud not be avoided here :-(
    integer(kind=4), allocatable :: cell_map(:)
    integer(kind=4), allocatable :: bnd_cell_map(:)

    ! Variables which follow are for backup saving to single file
    ! (These should probably be inside the Comm_Type)
    integer :: nc_s   ! number of cells in subdomain
    integer :: nb_s   ! number of bundary cells in subdoima
    integer :: nc_t   ! total number of cells 
    integer :: nb_t   ! total number of bundary cells

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
