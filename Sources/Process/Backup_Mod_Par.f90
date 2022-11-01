!==============================================================================!
  module Backup_Mod
!------------------------------------------------------------------------------!
!   Module containing functions to write and read backup files.                !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use User_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !-----------------!
  !   Backup type   !
  !-----------------!
  type Backup_Type
    integer :: interval
  end type
  type(Backup_Type) :: backup

  type(Mpi_File) :: fh

  contains

#   include "Backup_Mod/Load.f90"
#   include "Backup_Mod/Read_Cell_Real.f90"
#   include "Backup_Mod/Read_Int.f90"
#   include "Backup_Mod/Read_Int_Array.f90"
#   include "Backup_Mod/Read_Log.f90"
#   include "Backup_Mod/Read_Log_Array.f90"
#   include "Backup_Mod/Read_Real.f90"
#   include "Backup_Mod/Read_Real_Array.f90"
#   include "Backup_Mod/Read_Swarm.f90"
#   include "Backup_Mod/Read_Variable.f90"
#   include "Backup_Mod/Save.f90"
#   include "Backup_Mod/Time_To_Save.f90"
#   include "Backup_Mod/Write_Cell_Real.f90"
#   include "Backup_Mod/Write_Int.f90"
#   include "Backup_Mod/Write_Int_Array.f90"
#   include "Backup_Mod/Write_Log.f90"
#   include "Backup_Mod/Write_Log_Array.f90"
#   include "Backup_Mod/Write_Real.f90"
#   include "Backup_Mod/Write_Real_Array.f90"
#   include "Backup_Mod/Write_Swarm.f90"
#   include "Backup_Mod/Write_Variable.f90"

  end module
