!==============================================================================!
  module Info_Mod
!------------------------------------------------------------------------------!
!   Data types and functions for printing info boxes.                          !
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
  type System_Clock_Type
    integer :: cnt            ! system clock count rate
    integer :: ini, cur       ! system clock initial and current rate
    real    :: wall_time_max  ! maximum wall time
  end type

  !--------------------!
  !   Time_Info type   !
  !--------------------!
  type Time_Info_Type
    character(len=L_LINE) :: line_lead  = ''
    character(len=L_LINE) :: line_trail = ''
    character(len=L_LINE) :: lines(6)   = ''
  end type

  !--------------------!
  !   Iter_Info type   !
  !--------------------!
  type Iter_Info_Type
    integer               :: n_user_lines              = 0
    character(len=L_LINE) :: line_lead                 = ''
    character(len=L_LINE) :: line_sep                  = ''
    character(len=L_LINE) :: line_iter                 = ''
    character(len=L_LINE) :: lines(4)                  = ''
    character(len=L_LINE) :: lines_user(MAX_USER_LINES) = ''
  end type

  !--------------------!
  !   Bulk_Info type   !
  !--------------------!
  type Bulk_Info_Type
    character(len=L_LINE) :: line_lead  = ''
    character(len=L_LINE) :: line_foll  = ''
    character(len=L_LINE) :: line_sep   = ''
    character(len=L_LINE) :: line_trail = ''
    character(len=L_LINE) :: lines(3)   = ''
  end type

  type(System_Clock_Type), save :: sys_clock
  type(Time_Info_Type),    save :: time_info
  type(Iter_Info_Type),    save :: iter_info
  type(Bulk_Info_Type),    save :: bulk_info

  contains

#   include "Info_Mod/Start.f90"
#   include "Info_Mod/Time_To_Exit.f90"

#   include "Info_Mod/Time_Start.f90"
#   include "Info_Mod/Time_Fill.f90"
#   include "Info_Mod/Time_Print.f90"

#   include "Info_Mod/Iter_Start.f90"
#   include "Info_Mod/Iter_Fill.f90"
#   include "Info_Mod/Iter_Fill_At.f90"
#   include "Info_Mod/Iter_Fill_User_At.f90"
#   include "Info_Mod/Iter_Print.f90"

#   include "Info_Mod/Bulk_Start.f90"
#   include "Info_Mod/Bulk_Fill.f90"
#   include "Info_Mod/Bulk_Print.f90"

  end module
