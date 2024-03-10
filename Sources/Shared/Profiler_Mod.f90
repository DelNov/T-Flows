#include "../Shared/Assert.h90"

!==============================================================================!
  module Profiler_Mod
!------------------------------------------------------------------------------!
!>  The Profiler_Mod provides a framework for profiling the execution time of
!>  various functions or code segments in a program. It's designed to track
!>  and accumulate the time spent in individual functions, offering insights
!>  into the performance characteristics of the program.  It typically writes
!>  profiled results at the end of Convert and Process.
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Comm_Mod
  use Control_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  integer, parameter :: MAX_FUNCT = 2048  !! maximum number of functions or
                                          !! code segments profiler can track
  !-------------------!
  !   Profiler type   !
  !-------------------!
  !> Profiler_Type type encapsulates various
  !> components and procedures for profiling.
  type Profiler_Type

    integer, private :: n_functs                = 0  !! number of functions
                                                     !! being profiled
    integer, private :: curr_running            = 0  !! current function being
                                                     !! profiled.
    integer, private :: prev_running(MAX_FUNCT) = 0  !! storage for previously
                                                     !! running function indices

    character(DL), private :: funct_name(MAX_FUNCT)  !! stores the names of the
                                                     !! functions being profiled
    integer(DP),   private :: i_time_prev            !! system clock at prev
    integer(DP),   private :: i_time_curr            !! system clock at current
    integer(DP),   private :: sys_count_rate         !! system clock rate
    real,          private :: funct_time(MAX_FUNCT)  !! array to accumulate time
                                                     !!spent in each function
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
