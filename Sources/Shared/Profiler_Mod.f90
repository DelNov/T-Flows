!==============================================================================!
  module Profiler_Mod
!----------------------------------[Modules]-----------------------------------!
  use Comm_Mod
  use Control_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  integer, parameter :: MAX_FUNCT = 2048

  !-------------------!
  !   Profiler type   !
  !-------------------!
  type Profiler_Type

    integer, private :: n_functs                = 0
    integer, private :: curr_running            = 0
    integer, private :: prev_running(MAX_FUNCT) = 0

    character(DL), private :: funct_name(MAX_FUNCT)
    integer(DP),   private :: i_time_prev            ! system clock at prev
    integer(DP),   private :: i_time_curr            ! system clock at current
    integer(DP),   private :: sys_count_rate
    real,          private :: funct_time(MAX_FUNCT)  ! accumulated time

    contains
      procedure          :: Start
      procedure          :: Statistics
      procedure          :: Stop
      procedure, private :: Update_By_Rank

  end type

  !------------------------------------------!
  !   Create one instance of Profiler type   !
  !       for all other modules to use       !
  !------------------------------------------!
  type(Profiler_Type) :: Profiler

  contains

#   include "Profiler_Mod/Start.f90"
#   include "Profiler_Mod/Statistics.f90"
#   include "Profiler_Mod/Stop.f90"
#   include "Profiler_Mod/Update_By_Rank.f90"

  end module
