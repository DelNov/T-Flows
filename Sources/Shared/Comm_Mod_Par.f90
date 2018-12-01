!==============================================================================!
  module Comm_Mod
!------------------------------------------------------------------------------!
!   Module for MPI functionality.                                              !
!------------------------------------------------------------------------------!
  implicit none
!----------------------------------[Include]-----------------------------------!
  include 'mpif.h'
!==============================================================================!

  ! Parameters for size of typical variables in bytes
  integer, parameter :: SIZE_INT  = 8
  integer, parameter :: SIZE_REAL = 8

  integer :: this_proc  ! processor i.d.
  integer :: n_proc     ! number of processors

  ! Variables which follow are for backup saving to single file
  integer :: nc_s   ! number of cells in subdomain
  integer :: nb_s   ! number of bundary cells in subdoima
  integer :: nc_t   ! total number of cells 
  integer :: nb_t   ! total number of bundary cells

  integer :: cell_map_type
  integer :: bnd_cell_map_type

  contains

  include 'Comm_Mod/Parallel/Allocate.f90'
  include 'Comm_Mod/Parallel/Close_File.f90'
  include 'Comm_Mod/Parallel/Create_Buffers.f90'
  include 'Comm_Mod/Parallel/Create_New_Types.f90'
  include 'Comm_Mod/Parallel/End.f90'
  include 'Comm_Mod/Parallel/Exchange_Int.f90'
  include 'Comm_Mod/Parallel/Exchange_Real.f90'
  include 'Comm_Mod/Parallel/Global_Max_Real.f90'
  include 'Comm_Mod/Parallel/Global_Min_Real.f90'
  include 'Comm_Mod/Parallel/Global_Max_Int.f90'
  include 'Comm_Mod/Parallel/Global_Min_Int.f90'
  include 'Comm_Mod/Parallel/Global_Sum_Int_Array.f90'
  include 'Comm_Mod/Parallel/Global_Sum_Int.f90'
  include 'Comm_Mod/Parallel/Global_Sum_Real.f90'
  include 'Comm_Mod/Parallel/Load_Maps.f90'
  include 'Comm_Mod/Parallel/Open_File_Read.f90'
  include 'Comm_Mod/Parallel/Open_File_Write.f90'
  include 'Comm_Mod/Parallel/Read_Int.f90'
  include 'Comm_Mod/Parallel/Read_Bnd_Real.f90'
  include 'Comm_Mod/Parallel/Read_Cell_Real.f90'
  include 'Comm_Mod/Parallel/Read_Real.f90'
  include 'Comm_Mod/Parallel/Read_Text.f90'
  include 'Comm_Mod/Parallel/Set_Buffer_Color.f90'
  include 'Comm_Mod/Parallel/Start.f90'
  include 'Comm_Mod/Parallel/Wait.f90'
  include 'Comm_Mod/Parallel/Write_Int.f90'
  include 'Comm_Mod/Parallel/Write_Bnd_Real.f90'
  include 'Comm_Mod/Parallel/Write_Cell_Real.f90'
  include 'Comm_Mod/Parallel/Write_Real.f90'
  include 'Comm_Mod/Parallel/Write_Text.f90'

  end module
