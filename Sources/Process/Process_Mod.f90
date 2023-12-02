#include "../Shared/Assert.h90"
#include "../Shared/Browse.h90"
#include "../Shared/Unused.h90"

!==============================================================================!
  module Process_Mod
!------------------------------------------------------------------------------!
!>  The Process_Mod module is a centralized collection of functions crucial for
!>  conducting a CFD simulations. The module was primarily introduced to
!>  streamline the process of variable checking in Fortran and to enhance
!>  compatibility with the Intel compiler, especially in debug mode.  Still,
!>  it entails a wide range of procedures that are integral to various stages
!>  of CFD simulation. These include functions for initial variable set up,
!>  computing momentum (implicit and explicit), pressure, correcting velocity,
!>  computing energy, scalar variables, volume balancing, implementing the
!>  Rhie & Chow algorithm, PISO algorithm and updading boundary conditions.
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
  !> Encapsulation of procedures that manage the CFD solution algorithm.
  !> It is a framework which bundles together essential functions that
  !> manage various aspects of the CFD solution process.
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

  type(Process_Type) :: Process  !! singleton object of Process_Type, introduced
                                 !! to allow easy access to its procedures
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
