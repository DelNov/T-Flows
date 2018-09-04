!==============================================================================!
  module Monitor_Mod
!------------------------------------------------------------------------------!
!   Module for monitoring points.                                              !
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
 
  end type

  type(Monitor_Type), save :: monitor

  contains

  include 'Monitor_Mod/Initialize.f90'
  include 'Monitor_Mod/Finalize.f90'
  include 'Monitor_Mod/Write_4_Vars.f90'
  include 'Monitor_Mod/Write_5_Vars.f90'

  end module
