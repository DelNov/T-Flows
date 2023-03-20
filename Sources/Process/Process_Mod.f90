#include "../Shared/Assert.h90"
#include "../Shared/Browse.h90"
#include "../Shared/Unused.h90"

!==============================================================================!
  module Process_Mod
!------------------------------------------------------------------------------!
!   Collection of functions used in the Process program.  In honesty, it was   !
!   introduced to get rid of the Fortran header files with interfaces which,   !
!   in effect, was needed for Intel compiler to work in the debug mode.        !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Results_Mod
  use Read_Controls_Mod
  use Monitor_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !------------------!
  !   Process type   !
  !------------------!
  type Process_Type

    contains
      procedure :: Balance_Volume
      procedure :: Logo_Pro
      procedure :: Compute_Energy
      procedure :: Compute_Momentum_Explicit
      procedure :: Compute_Momentum
      procedure :: Compute_Pressure
      procedure :: Compute_Scalar
      procedure :: Convective_Outflow
      procedure :: Correct_Velocity
      procedure :: Initialize_Variables
      procedure :: Piso_Algorithm
      procedure :: Rhie_And_Chow
      procedure :: Update_Boundary_Values

  end type

  type(Process_Type) :: Process

  contains

#   include "Process_Mod/Balance_Volume.f90"
#   include "Process_Mod/Logo_Pro.f90"
#   include "Process_Mod/Compute_Energy.f90"
#   include "Process_Mod/Compute_Momentum_Explicit.f90"
#   include "Process_Mod/Compute_Momentum.f90"
#   include "Process_Mod/Compute_Pressure.f90"
#   include "Process_Mod/Compute_Scalar.f90"
#   include "Process_Mod/Convective_Outflow.f90"
#   include "Process_Mod/Correct_Velocity.f90"
#   include "Process_Mod/Initialize_Variables.f90"
#   include "Process_Mod/Piso_Algorithm.f90"
#   include "Process_Mod/Rhie_And_Chow.f90"
#   include "Process_Mod/Update_Boundary_Values.f90"

  end module
