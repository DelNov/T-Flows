!==============================================================================!
  module Backup_Mod
!------------------------------------------------------------------------------!
!   Module containing functions to write and read backup files.                !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
  use Comm_Mod
  use Name_Mod,    only: problem_name
  use Var_Mod
  use Turb_Mod
  use Field_Mod,   only: Field_Type, heat_transfer
  use Swarm_Mod,   only: Swarm_Mod_Find_Nearest_Node,  &
                         Swarm_Mod_Find_Nearest_Cell,  &
                         Swarm_Type, Particle_Type,    &
                         i_work, l_work, r_work, n_i_vars, n_l_vars, n_r_vars
  use Grid_Mod,    only: Grid_Type
  use Bulk_Mod,    only: Bulk_Type
  use Control_Mod
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
  include 'Backup_Mod/Read_Log_Array.f90'
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
  include 'Backup_Mod/Write_Log_Array.f90'
  include 'Backup_Mod/Write_Real.f90'
  include 'Backup_Mod/Write_Real_Array.f90'
  include 'Backup_Mod/Write_Swarm.f90'
  include 'Backup_Mod/Write_Variable.f90'

  end module
