# include "../Shared/Unused.h90"

!==============================================================================!
  module Generate_Mod
!------------------------------------------------------------------------------!
!>  Collection of functions used in the Generate program.  In full openness,
!>  it was introduced to get rid of the Fortran header files with interfaces
!>  which, in effect, was needed for Intel compiler to work in the debug mode.
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
  !> Contains (formerly global to Generate) routines for grid generation.
  type Generate_Type

    contains
      procedure :: Calculate_Geometry
      procedure :: Load_Dom
      procedure :: Logo_Gen
      procedure :: Print_Generate_Statistics
  end type

  type(Generate_Type) :: Generate  !! singleton type object Generate for easy
                                   !! (global) access to its member functions
  contains

#   include "Generate_Mod/Calculate_Geometry.f90"
#   include "Generate_Mod/Load_Dom.f90"
#   include "Generate_Mod/Logo_Gen.f90"
#   include "Generate_Mod/Print_Generate_Statistics.f90"

  end module
