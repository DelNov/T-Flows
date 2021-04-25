!==============================================================================!
  module Math_Mod
!------------------------------------------------------------------------------!
!   This module contains some basic mathematical (and related) functions       !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  ! Default precision for "Approx_Real" and "Smaller_Real" functions
  real, parameter :: DEFAULT_TOLERANCE = NANO

  contains

  include 'Math_Mod/Approx_Real.f90'
  include 'Math_Mod/Approx_String.f90'
  include 'Math_Mod/Cross_Product.f90'
  include 'Math_Mod/Distance.f90'
  include 'Math_Mod/Distance_Squared.f90'
  include 'Math_Mod/Gaussian_Elimination.f90'
  include 'Math_Mod/Harmonic_Mean.f90'
  include 'Math_Mod/Invert_3x3.f90'
  include 'Math_Mod/Rotate_Vector.f90'
  include 'Math_Mod/Smaller_Real.f90'
  include 'Math_Mod/Tet_Volume.f90'
  include 'Math_Mod/Tri_Area_Z.f90'

  end module
