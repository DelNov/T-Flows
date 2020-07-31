!==============================================================================!
  module Cpu_Timer_Mod
!----------------------------------[Modules]-----------------------------------!
  use Comm_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  integer, parameter :: MAX_FUNCT = 2048
  integer            :: n_funct   =    0
  integer            :: new_funct =    0
  integer            :: old_funct =    0

  character(SL) :: funct_name(MAX_FUNCT)
  real          :: funct_time(MAX_FUNCT)  ! accumulated time
  real          :: time_prev, time_curr

  contains

  include 'Cpu_Timer_Mod/Start.f90'
  include 'Cpu_Timer_Mod/Statistics.f90'
  include 'Cpu_Timer_Mod/Stop.f90'

  end module
