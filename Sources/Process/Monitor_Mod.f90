!==============================================================================!
  module Monitor_Mod
!------------------------------------------------------------------------------!
!   Module for monitoring points.                                              !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
  use Math_Mod
  use File_Mod
  use Comm_Mod
  use Grid_Mod
  use Field_Mod
  use Var_Mod
  use Control_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !------------------!
  !   Monitor type   !
  !------------------!
  type Monitor_Type

    type(Grid_Type), pointer :: pnt_grid  ! grid for which it is defined

    integer              :: n_points
    integer, allocatable :: cell(:)
    real,    allocatable :: x(:), y(:), z(:)
    integer, allocatable :: file_unit(:)

  end type

  contains

  include 'Monitor_Mod/Initialize.f90'
  include 'Monitor_Mod/Finalize.f90'
  include 'Monitor_Mod/Write_4_Vars.f90'
  include 'Monitor_Mod/Write_5_Vars.f90'

  end module
