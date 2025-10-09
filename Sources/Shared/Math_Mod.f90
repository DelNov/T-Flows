#include "../Shared/Unused.h90"

!==============================================================================!
  module Math_Mod
!------------------------------------------------------------------------------!
!>  This module is a holder of the Math_Type, which is a collection of
!>  some useful mathematical (and related) functions and subroutines used
!>  in T-Flows.
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  ! Default precision for "Approx_Real" and "Smaller_Real" functions
  real, parameter :: DEFAULT_TOLERANCE = NANO  !! default tolerance for
                                               !! Approx_Real and Smaller_Real
  !---------------!
  !   Math type   !
  !---------------!
  !> Math_Type contains a number of useful mathematical (and related) rountes.
  type Math_Type

    contains
      procedure          :: Approx_Real
      procedure          :: Approx_String
      procedure          :: Cross_Product
      procedure          :: Distance
      procedure          :: Distance_Squared
      procedure          :: Fit_Exp_Derivative_And_Two_Points
      procedure          :: Fit_Exp_Three_Points
      procedure          :: Gaussian_Elimination
      procedure          :: Harmonic_Mean
      procedure          :: Invert_3x3_Matrix
      procedure          :: Random_Real
      procedure          :: Rotate_Vector
      procedure, private :: Set_Array_Range
      procedure          :: Signed_Lower_Limit
      procedure          :: Signed_Upper_Limit
      procedure          :: Smaller_Real
      procedure          :: Tet_Inertia
      procedure          :: Tet_Volume
      procedure          :: Triangle_Area

  end type

  !--------------------------------------!
  !   Create one instance of Math type   !
  !     for all other modules to use     !
  !--------------------------------------!
  type(Math_Type) :: Math  !! definition of a singletone Math_Type,
                           !! in essence one global Math object

  contains

#   include "Math_Mod/Approx_Real.f90"
#   include "Math_Mod/Approx_String.f90"
#   include "Math_Mod/Cross_Product.f90"
#   include "Math_Mod/Distance.f90"
#   include "Math_Mod/Distance_Squared.f90"
#   include "Math_Mod/Fit_Exp_Derivative_And_Two_Points.f90"
#   include "Math_Mod/Fit_Exp_Three_Points.f90"
#   include "Math_Mod/Gaussian_Elimination.f90"
#   include "Math_Mod/Harmonic_Mean.f90"
#   include "Math_Mod/Invert_3x3_Matrix.f90"
#   include "Math_Mod/Random_Real.f90"
#   include "Math_Mod/Rotate_Vector.f90"
#   include "Math_Mod/Set_Array_Range.f90"
#   include "Math_Mod/Signed_Lower_Limit.f90"
#   include "Math_Mod/Signed_Upper_Limit.f90"
#   include "Math_Mod/Smaller_Real.f90"
#   include "Math_Mod/Tet_Inertia.f90"
#   include "Math_Mod/Tet_Volume.f90"
#   include "Math_Mod/Triangle_Area.f90"

  end module
