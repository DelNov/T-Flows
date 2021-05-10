!==============================================================================!
  module Cpu_Timer_Mod
!----------------------------------[Modules]-----------------------------------!
  use Comm_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  integer, parameter :: MAX_FUNCT = 2048

  !--------------------!
  !   Cpu_Timer type   !
  !--------------------!
  type Cpu_Timer_Type

    integer, private :: n_funct   =    0
    integer, private :: new_funct =    0
    integer, private :: old_funct =    0

    character(DL), private :: funct_name(MAX_FUNCT)
    real,          private :: funct_time(MAX_FUNCT)  ! accumulated time
    real,          private :: time_prev, time_curr

    contains
      procedure :: Start
      procedure :: Statistics
      procedure :: Stop

  end type

  !-------------------------------------------!
  !   Create one instance of Cpu_Timer type   !
  !       for all other modules to use        !
  !-------------------------------------------!
  type(Cpu_Timer_Type) :: Cpu_Timer

  contains

  include 'Cpu_Timer_Mod/Start.f90'
  include 'Cpu_Timer_Mod/Statistics.f90'
  include 'Cpu_Timer_Mod/Stop.f90'

  end module
