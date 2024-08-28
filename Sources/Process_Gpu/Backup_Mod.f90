#include "../Shared/Assert.h90"
#include "../Shared/Browse.h90"
#include "../Shared/Unused.h90"

!==============================================================================!
  module Backup_Mod
!------------------------------------------------------------------------------!
!   Module containing functions to write and read backup files.                !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Time_Mod
  use Turb_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !-----------------!
  !   Backup type   !
  !-----------------!
  type Backup_Type
    integer :: interval

    contains
      procedure :: Load
      procedure :: Load_Cell_Real
      procedure :: Load_Int
      procedure :: Load_Int_Array
      procedure :: Load_Log
      procedure :: Load_Log_Array
      procedure :: Load_Real
      procedure :: Load_Real_Array
      procedure :: Load_Variable
      procedure :: Save
      procedure :: Save_Cell_Real
      procedure :: Save_Int
      procedure :: Save_Int_Array
      procedure :: Save_Log
      procedure :: Save_Log_Array
      procedure :: Save_Real
      procedure :: Save_Real_Array
      procedure :: Save_Variable
      procedure :: Time_To_Save_Backup

  end type

  !----------------------!
  !   Singleton object   !
  !----------------------!
  type(Backup_Type) :: Backup

# if T_FLOWS_MPI == 1
  type(Mpi_File) :: fh
# else
  integer :: fh
# endif

  contains

#   include "Backup_Mod/Load.f90"
#   include "Backup_Mod/Load_Cell_Real.f90"
#   include "Backup_Mod/Load_Int.f90"
#   include "Backup_Mod/Load_Int_Array.f90"
#   include "Backup_Mod/Load_Log.f90"
#   include "Backup_Mod/Load_Log_Array.f90"
#   include "Backup_Mod/Load_Real.f90"
#   include "Backup_Mod/Load_Real_Array.f90"
#   include "Backup_Mod/Load_Variable.f90"
#   include "Backup_Mod/Save.f90"
#   include "Backup_Mod/Save_Cell_Real.f90"
#   include "Backup_Mod/Save_Int.f90"
#   include "Backup_Mod/Save_Int_Array.f90"
#   include "Backup_Mod/Save_Log.f90"
#   include "Backup_Mod/Save_Log_Array.f90"
#   include "Backup_Mod/Save_Real.f90"
#   include "Backup_Mod/Save_Real_Array.f90"
#   include "Backup_Mod/Save_Variable.f90"
#   include "Backup_Mod/Time_To_Save_Backup.f90"

  end module
