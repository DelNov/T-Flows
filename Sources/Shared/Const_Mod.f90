!==============================================================================!
  module Const_Mod
!------------------------------------------------------------------------------!
!>  Centralized repository for various constants used across different parts
!>  of T-Flows.  It defines a wide range of constants including program names,
!>  standard string lengths, precision types, mathematical constants, and
!>  parameters for multi-domain simulations  It is used by all other modules.
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  ! Program name (T_FLOWS_PROGRAM) is passed from makefile of each program
# if T_FLOWS_PROGRAM == 1
  character(7), parameter :: PROGRAM_NAME = "Convert"
# elif T_FLOWS_PROGRAM == 2
  character(8), parameter :: PROGRAM_NAME = "Generate"
# elif T_FLOWS_PROGRAM == 3
  character(6), parameter :: PROGRAM_NAME = "Divide"
# elif T_FLOWS_PROGRAM == 4
  character(7), parameter :: PROGRAM_NAME = "Process"
# else
  character(9), parameter :: PROGRAM_NAME = "Undefined"
# endif

  ! Standard string length
  integer, parameter :: VL =   4
    !! name length for T-Flows variables (such as 'u', 'v', 'w', 'kin, 'vof',
    !! 'c_01, 'q_01', ...) used to identify them when specifying boundary and
    !! initial conditions in the control file, and when reporting their
    !! convergence characteristics on the terminal through Info_Mod
  integer, parameter :: SL =  80  !! standard string length (like page width)
  integer, parameter :: DL = 160  !! double string length (twice the page width)

  ! Maximum number of string items in a line
  integer, parameter :: MAX_STRING_ITEMS = 32
    !! maximum number of string items in a line used mostly to read options
    !! for PETSc solvers, boundary and initial conditions, or porosity regions

  ! Double and single precision constants definitions
  integer, parameter :: DP = 8  !! double precisions for real and long integer
  integer, parameter :: SP = 4  !! single precisions for real and short integer

  ! Precision of integers and floating point numbers.  These parameters
  ! assume values according to what is in the makefile (which is a bad
  ! practice but that's another story).  But that's not enough for program
  ! to work.  All MPI-based routines would need to distinguish MPI_INTEGER
  ! from MPI_INTEGER8.  Also, you will have to link the program with proper
  ! METIS library, check all the makefiles for details on that.
  integer, parameter :: IP = sizeof(1)   !! integer precision
  integer, parameter :: LP = IP          !! logical precision (always the same
                                         !! as integer precision in Fortran)
  integer, parameter :: RP = sizeof(1.0) !! real number precision (in this form
                                         !! it is inherited from makefiles)

  ! Version of the .cfn and .dim files
  ! (First four digits are the year, last two the month)
  integer, parameter :: VERSION_CFN    = 202310  !! version of the .cfn files
  integer, parameter :: VERSION_DIM    = 202304  !! version of the .dim files
  integer, parameter :: VERSION_BACKUP = 202304  !! version of the backup files

  !----------------------------------------!
  !   A few handy mathematical constants   !
  !----------------------------------------!

  ! Big and small numbers in metric system to avoid ghost numbers
  real, parameter :: YOTTA = 1.e+24  !! avoid ghost number 1.0e+24
  real, parameter :: ZETTA = 1.e+21  !! avoid ghost number 1.0e+21
  real, parameter :: EXA   = 1.e+18  !! avoid ghost number 1.0e+18
  real, parameter :: PETA  = 1.e+15  !! avoid ghost number 1.0e+15
  real, parameter :: TERA  = 1.e+12  !! avoid ghost number 1.0e+12
  real, parameter :: GIGA  = 1.e+9   !! avoid ghost number 1.0e+9
  real, parameter :: MEGA  = 1.e+6   !! avoid ghost number 1.0e+6
  real, parameter :: KILO  = 1.e+3   !! avoid ghost number 1.0e+3
  real, parameter :: MILI  = 1.e-3   !! avoid ghost number 1.0e-3
  real, parameter :: MICRO = 1.e-6   !! avoid ghost number 1.0e-6
  real, parameter :: NANO  = 1.e-9   !! avoid ghost number 1.0e-9
  real, parameter :: PICO  = 1.e-12  !! avoid ghost number 1.0e-12
  real, parameter :: FEMTO = 1.e-15  !! avoid ghost number 1.0e-15
  real, parameter :: ATTO  = 1.e-18  !! avoid ghost number 1.0e-18
  real, parameter :: ZEPTO = 1.e-21  !! avoid ghost number 1.0e-21
  real, parameter :: YOCTO = 1.e-24  !! avoid ghost number 1.0e-24

  real,    parameter :: HUGE     = PETA        !! a very big (huge) number
  real,    parameter :: TINY     = FEMTO       !! a very small (tiny) number
  integer, parameter :: HUGE_INT = 1073741824  !! big integer (this is 2^30)

  ! Euler's prime number (also the largest integer in 32 bit precision)
  integer, parameter :: EULER    = 2147483647  ! Euler's prime number 2^31 - 1

  ! Archimedes’ constant
  real, parameter :: PI = 3.14159265359  !! Archimedes constant

  ! These are often used in turbulence models
  real, parameter :: ONE_THIRD  = 1.0 / 3.0        !! avoids frequent 1.0/3.0
  real, parameter :: TWO_THIRDS = 1.0 - ONE_THIRD  !! avoids frequent 2.0/3.0
  real, parameter :: ONE_SIXTH  = ONE_THIRD * 0.5  !! avoids frequent 1.0/6.0

  !------------------------------------------------------!
  !   Constants related to multiple domain simulations   !
  !------------------------------------------------------!

  ! Maximum number of domains
  integer, parameter :: MD = 4  !! maximum number of problem domains

  ! Maximum number of variables exchanged at the interface
  integer, parameter :: MAX_VARS_INTERFACE = 9  !! maximum number of variables
                                                !! exchanged at the interface
                                                !! between two domains
  end module
