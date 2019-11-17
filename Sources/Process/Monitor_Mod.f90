!==============================================================================!
  module Monitor_Mod
!------------------------------------------------------------------------------!
!   Module for monitoring points.                                              !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod,   only: HUGE, TINY
  use Math_Mod
  use File_Mod
  use Comm_Mod,    only: n_proc, Comm_Mod_Global_Min_Real
  use Grid_Mod,    only: Grid_Type
  use Field_Mod,   only: Field_Type
  use Var_Mod,     only: Var_Type
  use Control_Mod, only: Control_Mod_Read_Int_Item,  &
                         Control_Mod_Read_Real_Array
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !------------------!
  !   Monitor type   !
  !------------------!
  type Monitor_Type

    integer              :: n_points
    integer, allocatable :: cell(:)
    real,    allocatable :: x(:), y(:), z(:)
    integer, allocatable :: file_unit(:)

  end type

  type(Monitor_Type), save :: monitor

  contains

  include 'Monitor_Mod/Initialize.f90'
  include 'Monitor_Mod/Finalize.f90'
  include 'Monitor_Mod/Write_4_Vars.f90'
  include 'Monitor_Mod/Write_5_Vars.f90'

  end module
