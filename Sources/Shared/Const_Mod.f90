!==============================================================================!
  module Const_Mod
!------------------------------------------------------------------------------!
!   Constants definitions for all other modules.                               !
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !----------------------------------------!
  !   A few handy mathematical constants   !
  !----------------------------------------!
  real,    parameter :: HUGE     = 1.e+30
  real,    parameter :: TINY     = 1.e-30
  integer, parameter :: HUGE_INT = 9223372036854775807

  ! Big and small numbers in metric system
  real, parameter :: TERA  = 1.e+12
  real, parameter :: GIGA  = 1.e+9
  real, parameter :: MEGA  = 1.e+6
  real, parameter :: KILO  = 1.e+3
  real, parameter :: MILI  = 1.e-3
  real, parameter :: MICRO = 1.e-6
  real, parameter :: NANO  = 1.e-9
  real, parameter :: PICO  = 1.e-12

  ! Archimedes’ constant
  real, parameter :: PI = 3.14159265359

  ! Gravitational constant for Earth (where we are most of the time)
  real, parameter :: EARTH_G = 9.807

  ! These are often used in turbulence models
  real, parameter :: ONE_THIRD  = 0.33333333333333333
  real, parameter :: TWO_THIRDS = 1.0 - ONE_THIRD
  real, parameter :: ONE_SIXTH  = ONE_THIRD * 0.5

  end module 
