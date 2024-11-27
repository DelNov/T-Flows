!==============================================================================!
  module Info_Mod
!------------------------------------------------------------------------------!
!>  The Info_Mod module in T-Flows is designed for creating formatted output
!>  related to simulation progress.  It appears organizes and display various
!>  types of information, including system clock data, time step details,
!>  iteration updates, residuals achieved in solution of linear system of
!>  equation as well as the outer (SIMPLE/PISO) iteration residuals.
!>  In addition, it also prints bulk flow information including volume flow
!>  rates and pressure drops in each coordiate direction.
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Field_Mod
!------------------------------------------------------------------------------!
  implicit none
!-----------------------------------[Locals]-----------------------------------!
  integer, parameter :: L_LINE         = 127
  integer, parameter :: L_BOX          =  21
  integer, parameter :: MAX_USER_LINES = 256
!==============================================================================!

  !-----------------------!
  !   System Clock type   !
  !-----------------------!
  !> Manages the system clock count rate and wall time, useful
  !> for tracking real-time performance of the simulation.
  type System_Clock_Type
    integer(DP) :: cnt            !! system clock count rate
    integer(DP) :: ini, cur       !! system clock initial and current rate
    real        :: wall_time_max  !! maximum wall time
  end type

  !--------------------!
  !   Time Info type   !
  !--------------------!
  !> Responsible for formatting and displaying information about the
  !> simulation's time steps. This includes current time step, physical
  !> time, and wall-clock time.
  type Time_Info_Type
    character(len=L_LINE) :: line_lead  = ''
    character(len=L_LINE) :: line_trail = ''
    character(len=L_LINE) :: line(6)    = ''
  end type

  !--------------------!
  !   Iter Info type   !
  !--------------------!
  !> Manages the display of information related to iterations within each
  !> time step. It includes functionality to handle customized user lines,
  !> allowing for extended information display about the iterative process.
  type Iter_Info_Type
    integer               :: n_user_lines               = 0
    character(len=L_LINE) :: line_lead                  = ''
    character(len=L_LINE) :: line_sep                   = ''
    character(len=L_LINE) :: line_iter                  = ''
    character(len=L_LINE) :: line(4)                    = ''
    character(len=L_LINE) :: lines_user(MAX_USER_LINES) = ''
  end type

  !--------------------!
  !   Bulk Info type   !
  !--------------------!
  !> Handles the display of bulk flow information such as volume flow rates
  !> and pressure drops in different directions. This is crucial for
  !> understanding the overall flow behavior in the simulation domain.
  type Bulk_Info_Type
    character(len=L_LINE) :: line_lead  = ''
    character(len=L_LINE) :: line_foll  = ''
    character(len=L_LINE) :: line_sep   = ''
    character(len=L_LINE) :: line_trail = ''
    character(len=L_LINE) :: line(3)    = ''
  end type

  !---------------!
  !   Info type   !
  !---------------!
  !> A composite type that encapsulates the system clock, time, iteration,
  !> and bulk information types. It contains procedures for starting,
  !> filling, and printing information related to time steps, iterations,
  !> and bulk flow parameters.
  type Info_Type
    type(System_Clock_Type), private :: clock
    type(Time_Info_Type),    private :: time
    type(Iter_Info_Type),    private :: iter
    type(Bulk_Info_Type),    private :: bulk

    contains
      procedure :: Start_Info
      procedure :: Time_To_Exit

      procedure :: Time_Start
      procedure :: Time_Fill
      procedure :: Time_Print

      procedure :: Iter_Start
      procedure :: Iter_Fill
      procedure :: Iter_Fill_At
      procedure :: Iter_Fill_Scalar_At
      procedure :: Iter_Print

      procedure :: Bulk_Start
      procedure :: Bulk_Fill
      procedure :: Bulk_Print

  end type

  !---------------------------------!
  !   Singletone type Info object   !
  !---------------------------------!
  !> A globally accessible instance of the Info_Type, allowing for convenient
  !> access and manipulation of simulation information throughout T-Flows.
  type(Info_Type) :: Info

  contains

#   include "Info_Mod/Start_Info.f90"
#   include "Info_Mod/Time_To_Exit.f90"

#   include "Info_Mod/Time_Start.f90"
#   include "Info_Mod/Time_Fill.f90"
#   include "Info_Mod/Time_Print.f90"

#   include "Info_Mod/Iter_Start.f90"
#   include "Info_Mod/Iter_Fill.f90"
#   include "Info_Mod/Iter_Fill_At.f90"
#   include "Info_Mod/Iter_Fill_Scalar_At.f90"
#   include "Info_Mod/Iter_Print.f90"

#   include "Info_Mod/Bulk_Start.f90"
#   include "Info_Mod/Bulk_Fill.f90"
#   include "Info_Mod/Bulk_Print.f90"

  end module
