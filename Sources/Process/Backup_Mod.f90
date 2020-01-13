!==============================================================================!
  module Backup_Mod
!------------------------------------------------------------------------------!
!   Module containing functions to write and read backup files.                !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
  use Cpu_Timer_Mod
  use Comm_Mod
  use File_Mod
  use Var_Mod
  use Turb_Mod
  use Field_Mod
  use Swarm_Mod
  use Multiphase_Mod
  use Grid_Mod
  use Bulk_Mod
  use Control_Mod
  use User_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  contains

  include 'Backup_Mod/Load.f90'
  include 'Backup_Mod/Read_Bnd.f90'
  include 'Backup_Mod/Read_Cell_Bnd.f90'
  include 'Backup_Mod/Read_Cell.f90'
  include 'Backup_Mod/Read_Face.f90'
  include 'Backup_Mod/Read_Int.f90'
  include 'Backup_Mod/Read_Int_Array.f90'
  include 'Backup_Mod/Read_Real.f90'
  include 'Backup_Mod/Read_Real_Array.f90'
  include 'Backup_Mod/Read_Swarm.f90'
  include 'Backup_Mod/Read_Variable.f90'
  include 'Backup_Mod/Save.f90'
  include 'Backup_Mod/Write_Bnd.f90'
  include 'Backup_Mod/Write_Cell_Bnd.f90'
  include 'Backup_Mod/Write_Cell.f90'
  include 'Backup_Mod/Write_Face.f90'
  include 'Backup_Mod/Write_Int.f90'
  include 'Backup_Mod/Write_Int_Array.f90'
  include 'Backup_Mod/Write_Real.f90'
  include 'Backup_Mod/Write_Real_Array.f90'
  include 'Backup_Mod/Write_Swarm.f90'
  include 'Backup_Mod/Write_Variable.f90'

  end module
