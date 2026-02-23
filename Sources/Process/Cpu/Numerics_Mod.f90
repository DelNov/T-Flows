#include "../../Shared/Browse.h90"

!==============================================================================!
  module Numerics_Mod
!------------------------------------------------------------------------------!
!>  The Numerics_Mod module in T-Flows is essential for implementing numerical
!>  strategies used in the discretization of conservation equations via the
!>  finite volume method. It centralizes parameters and procedures common to
!>  various discretization schemes, reducing code duplication and improving
!>  maintainability. This module does not handle linear solvers, which are
!>  managed separately in modules like Native_Mod and Petsc_Mod. Numerics_Mod
!>  includes advection schemes, time integration parameters, algorithms for
!>  pressure-velocity coupling, and methods for gradient computation.  It does
!>  not define a new type.  Rather, it is a collection of procedures (whose
!>  names start with Numerics_Mod) and parameters used when dicretizing
!>  the equations.
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Matrix_Mod
  use Var_Mod
  use Work_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  ! Parameters for advection scheme
  integer, parameter :: UPWIND    = 40009
  integer, parameter :: CENTRAL   = 40013
  integer, parameter :: LUDS      = 40031
  integer, parameter :: QUICK     = 40037
  integer, parameter :: SMART     = 40039
  integer, parameter :: GAMMA     = 40063
  integer, parameter :: MINMOD    = 40087
  integer, parameter :: BLENDED   = 40093
  integer, parameter :: SUPERBEE  = 40099
  integer, parameter :: AVL_SMART = 40111
  integer, parameter :: CICSAM    = 40123
  integer, parameter :: STACS     = 40127

  ! Time integration parameters
  integer, parameter :: LINEAR        = 40129
  integer, parameter :: PARABOLIC     = 40151

  ! Algorithms for pressure velocity coupling
  integer, parameter :: SIMPLE = 40153
  integer, parameter :: PISO   = 40163
  integer, parameter :: CHOI   = 40169

  ! Gradient computation algorithms
  integer, parameter :: LEAST_SQUARES = 40177
  integer, parameter :: GAUSS_THEOREM = 40189

  contains

#   include "Numerics_Mod/Advection_Scheme.f90"
#   include "Numerics_Mod/Advection_Scheme_Code.f90"
#   include "Numerics_Mod/Advection_Term.f90"
#   include "Numerics_Mod/Gradient_Method_Code.f90"
#   include "Numerics_Mod/Inertial_Term.f90"
#   include "Numerics_Mod/Min_Max.f90"
#   include "Numerics_Mod/Pressure_Momentum_Coupling_Code.f90"
#   include "Numerics_Mod/Time_Integration_Scheme_Code.f90"
#   include "Numerics_Mod/Under_Relax.f90"

  end module
