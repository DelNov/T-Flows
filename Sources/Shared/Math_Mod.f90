!==============================================================================!
  module Math_Mod
!------------------------------------------------------------------------------!
!   This is a prototype of a small module which would contain some basic       !
!   mathematic and related functions.                                          !
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  ! Default precision for "Approx_Real" and "Smaller_Real" functions
  real, parameter :: DEFAULT_TOLERANCE = 1.0e-9

  contains

  include 'Math_Mod/Approx_Real.f90'
  include 'Math_Mod/Approx_String.f90'
  include 'Math_Mod/Cross_Product.f90'
  include 'Math_Mod/Distance.f90'
  include 'Math_Mod/Distance_Squared.f90'
  include 'Math_Mod/Smaller_Real.f90'
  include 'Math_Mod/Tet_Volume.f90'
  include 'Math_Mod/Tri_Area_Z.f90'

  end module
