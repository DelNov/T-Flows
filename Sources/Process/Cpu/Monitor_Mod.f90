#include "../../Shared/Browse.h90"

!==============================================================================!
  module Monitor_Mod
!------------------------------------------------------------------------------!
!>  This module, Monitor_Mod, is dedicated to managing monitoring points in
!>  the Process sub-program of T-Flows. It enables the tracking and reporting
!>  of selected variables (velocity components, pressure, temperature, and
!>  scalars) at specified points within the computational domain. The
!>  monitoring points are user-defined and the module handles their data
!>  storage, initialization, finalization, and writing the monitored variable
!>  data to specific monitoring files.
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Field_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !------------------!
  !   Monitor type   !
  !------------------!
  type Monitor_Type

    type(Grid_Type), pointer :: pnt_grid  !! grid for which it is defined

    integer              :: n_points          !! number of momnitoring points
    integer, allocatable :: cell(:)           !! nearest cell to each point
    real,    allocatable :: x(:), y(:), z(:)  !! points' coordinate
    integer, allocatable :: file_unit(:)      !! unit for monitoring file

    contains
      procedure :: Initialize
      procedure :: Finalize
      procedure :: Write_Vars

  end type

  contains

#   include "Monitor_Mod/Initialize.f90"
#   include "Monitor_Mod/Finalize.f90"
#   include "Monitor_Mod/Write_Vars.f90"

  end module
