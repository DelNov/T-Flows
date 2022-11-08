!==============================================================================!
  module Generate_Mod
!------------------------------------------------------------------------------!
!   Collection of functions used in the Generate program.  In honesty, it was   !
!   introduced to get rid of the Fortran header files with interfaces which,   !
!   in effect, was needed for Intel compiler to work in the debug mode.        !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Domain_Mod
  use Smooths_Mod
  use Refines_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !-------------------!
  !   Generate type   !
  !-------------------!
  type Generate_Type

    contains
      procedure :: Calculate_Geometry
      procedure :: Load_Dom
      procedure :: Logo_Gen
  end type

  type(Generate_Type) :: Generate

  contains

#   include "Generate_Mod/Calculate_Geometry.f90"
#   include "Generate_Mod/Load_Dom.f90"
#   include "Generate_Mod/Logo_Gen.f90"

  end module
