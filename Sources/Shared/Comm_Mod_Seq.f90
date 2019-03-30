!==============================================================================!
  module Comm_Mod
!------------------------------------------------------------------------------!
!   Module for no MPI functionality.                                           !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod, only: Grid_Type
!------------------------------------------------------------------------------!
  implicit none
!----------------------------------[Include]-----------------------------------!
! include 'mpif.h'
!==============================================================================!

  ! Parameters for size of typical variables in bytes
  integer, parameter :: SIZE_INT  = 8
  integer, parameter :: SIZE_LOG  = 8
  integer, parameter :: SIZE_REAL = 8

  integer :: this_proc  ! processor i.d.
  integer :: n_proc     ! number of processors

  ! Variables which follow are for backup saving to single file
  integer :: nc_s   ! number of cells in subdomain
  integer :: nb_s   ! number of bundary cells in subdoima
  integer :: nc_t   ! total number of cells 
  integer :: nb_t   ! total number of bundary cells

  contains

  include 'Comm_Mod/Sequential/Allocate.f90'
  include 'Comm_Mod/Sequential/Close_File.f90'
  include 'Comm_Mod/Sequential/Create_Buffers.f90'
  include 'Comm_Mod/Sequential/Create_New_Types.f90'
  include 'Comm_Mod/Sequential/End.f90'
  include 'Comm_Mod/Sequential/Exchange_Int.f90'
  include 'Comm_Mod/Sequential/Exchange_Real.f90'
  include 'Comm_Mod/Sequential/Global_Max_Real.f90'
  include 'Comm_Mod/Sequential/Global_Min_Real.f90'
  include 'Comm_Mod/Sequential/Global_Max_Int.f90'
  include 'Comm_Mod/Sequential/Global_Min_Int.f90'
  include 'Comm_Mod/Sequential/Global_Sum_Int.f90'
  include 'Comm_Mod/Sequential/Global_Sum_Int_Array.f90'
  include 'Comm_Mod/Sequential/Global_Sum_Real.f90'
  include 'Comm_Mod/Sequential/Global_Sum_Real_Array.f90'
  include 'Comm_Mod/Sequential/Load_Maps.f90'
  include 'Comm_Mod/Sequential/Open_File_Read.f90'
  include 'Comm_Mod/Sequential/Open_File_Write.f90'
  include 'Comm_Mod/Sequential/Read_Int.f90'
  include 'Comm_Mod/Sequential/Read_Int_Array.f90'
  include 'Comm_Mod/Sequential/Read_Log_Array.f90'
  include 'Comm_Mod/Sequential/Read_Bnd_Real.f90'
  include 'Comm_Mod/Sequential/Read_Cell_Real.f90'
  include 'Comm_Mod/Sequential/Read_Real.f90'
  include 'Comm_Mod/Sequential/Read_Real_Array.f90'
  include 'Comm_Mod/Sequential/Read_Text.f90'
  include 'Comm_Mod/Sequential/Set_Buffer_Color.f90'
  include 'Comm_Mod/Sequential/Start.f90'
  include 'Comm_Mod/Sequential/Wait.f90'
  include 'Comm_Mod/Sequential/Write_Int.f90'
  include 'Comm_Mod/Sequential/Write_Int_Array.f90'
  include 'Comm_Mod/Sequential/Write_Log_Array.f90'
  include 'Comm_Mod/Sequential/Write_Bnd_Real.f90'
  include 'Comm_Mod/Sequential/Write_Cell_Real.f90'
  include 'Comm_Mod/Sequential/Write_Real.f90'
  include 'Comm_Mod/Sequential/Write_Real_Array.f90'
  include 'Comm_Mod/Sequential/Write_Text.f90'

  end module
