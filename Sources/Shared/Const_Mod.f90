!==============================================================================!
  module Const_Mod
!------------------------------------------------------------------------------!
!   Constants definitions for all other modules.                               !
!                                                                              !
!   Program name (T_FLOWS_PROGRAM) is passed from makefile of each program     !
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

# if T_FLOWS_PROGRAM == 1
  character(8), parameter :: PROGRAM_NAME = "Generate"
# elif T_FLOWS_PROGRAM == 2
  character(6), parameter :: PROGRAM_NAME = "Divide"
# elif T_FLOWS_PROGRAM == 3
  character(7), parameter :: PROGRAM_NAME = "Convert"
# elif T_FLOWS_PROGRAM == 4
  character(7), parameter :: PROGRAM_NAME = "Process"
# else
  character(9), parameter :: PROGRAM_NAME = "Undefined"
# endif

  ! Standard string length
  integer, parameter :: VL =   4  ! variable name length
  integer, parameter :: SL =  80  ! standard string length
  integer, parameter :: DL = 160  ! double string length

  ! Maximum number of string items in a line
  ! (for example when reading boundary conditions, options for PETSc)
  integer, parameter :: MAX_STRING_ITEMS = 32

  ! Double and single precision constants definitions
  integer, parameter :: DP =  8  ! double precisions for real and integer
  integer, parameter :: SP =  4  ! single precisions for real and integer

  ! Precision of integers and floating point numbers.  These parameters
  ! assume values according to what is in the makefile (which is a bad
  ! practice but that's another story).  But that's not enough for program
  ! to work.  All MPI-based routines would need to have MPI_INTEGER8
  ! instead of MPI_INTEGER.
  ! Also, you will have to link the program with proper METIS library, check
  ! all the makefiles for details on that.
  integer, parameter :: IP = sizeof(1)   ! integer precision
  integer, parameter :: LP = IP          ! logical precision
  integer, parameter :: RP = sizeof(1.0) ! real number precision

  ! Version of the .cfn and .dim files
  ! (This number was set in April of 2023)
  integer, parameter :: VERSION_CFN    = 202304
  integer, parameter :: VERSION_DIM    = 202304
  integer, parameter :: VERSION_BACKUP = 202304

  !----------------------------------------!
  !   A few handy mathematical constants   !
  !----------------------------------------!

  ! Big and small numbers in metric system
  real, parameter :: YOTTA = 1.e+24
  real, parameter :: ZETTA = 1.e+21
  real, parameter :: EXA   = 1.e+18
  real, parameter :: PETA  = 1.e+15
  real, parameter :: TERA  = 1.e+12
  real, parameter :: GIGA  = 1.e+9
  real, parameter :: MEGA  = 1.e+6
  real, parameter :: KILO  = 1.e+3
  real, parameter :: MILI  = 1.e-3
  real, parameter :: MICRO = 1.e-6
  real, parameter :: NANO  = 1.e-9
  real, parameter :: PICO  = 1.e-12
  real, parameter :: FEMTO = 1.e-15
  real, parameter :: ATTO  = 1.e-18
  real, parameter :: ZEPTO = 1.e-21
  real, parameter :: YOCTO = 1.e-24

  real,    parameter :: HUGE     = PETA
  real,    parameter :: TINY     = FEMTO
  integer, parameter :: HUGE_INT = 1073741824  ! this is 2 ** 30

  ! Euler's prime number (also the largest integer in 32 bit precision)
  integer, parameter :: EULER    = 2147483647  ! 2 ** 31 - 1

  ! Archimedes’ constant
  real, parameter :: PI = 3.14159265359

  ! These are often used in turbulence models
  real, parameter :: ONE_THIRD  = 1.0 / 3.0
  real, parameter :: TWO_THIRDS = 1.0 - ONE_THIRD
  real, parameter :: ONE_SIXTH  = ONE_THIRD * 0.5

  !------------------------------------------------------!
  !   Constants related to multiple domain simulations   !
  !------------------------------------------------------!

  ! Maximum number of domains
  integer, parameter :: MD  = 4

  ! Maximum number of variables exchanged at the interface
  integer, parameter :: MAX_VARS_INTERFACE = 9

  end module
