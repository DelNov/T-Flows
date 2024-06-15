#include "../Shared/Assert.h90"
#include "../Shared/Browse.h90"
#include "../Shared/Unused.h90"

! A couple of macros to shorten the syntax a bit
#define Inc(X, Y)  X = X + Y
#define Dec(X, Y)  X = X - Y
#define Div(X, Y)  X = X / Y
#define Mul(X, Y)  X = X * Y

!==============================================================================!
  module Process_Mod
!------------------------------------------------------------------------------!
  use Assert_Mod
  use Profiler_Mod
  use Iter_Mod
  use Gpu_Mod
  use Control_Mod
  use Results_Mod
  use Read_Controls_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !------------------!
  !                  !
  !   Process type   !
  !                  !
  !------------------!
  type Process_Type

    contains

      ! General conservation equation
      procedure :: Form_System_Matrix
      procedure :: Add_Advection_Term
      procedure :: Add_Inertial_Term

      ! Related to momentum conservation
      procedure :: Compute_Momentum
      procedure :: Insert_Momentum_Bc
      procedure :: Add_Pressure_Term
      procedure :: Correct_Velocity
      procedure :: Initialize_Variables
      procedure :: Update_Boundary_Values

      ! Related to pressure solution
      procedure :: Compute_Pressure
      procedure :: Form_Pressure_Matrix
      procedure :: Insert_Volume_Source_For_Pressure

      ! Energy / Enthalpy
      procedure :: Compute_Energy
      procedure :: Insert_Energy_Bc

  end type

  type(Process_Type) :: Process  !! singleton object of Process_Type, introduced
                                 !! to allow easy access to its procedures
  contains

    ! General conservation equation
#   include "Process_Mod/Form_System_Matrix.f90"
#   include "Process_Mod/Add_Advection_Term.f90"
#   include "Process_Mod/Add_Inertial_Term.f90"

    ! Related to momentum conservation
#   include "Process_Mod/Compute_Momentum.f90"
#   include "Process_Mod/Insert_Momentum_Bc.f90"
#   include "Process_Mod/Add_Pressure_Term.f90"
#   include "Process_Mod/Correct_Velocity.f90"
#   include "Process_Mod/Initialize_Variables.f90"
#   include "Process_Mod/Update_Boundary_Values.f90"

    ! Related to pressure solution
#   include "Process_Mod/Compute_Pressure.f90"
#   include "Process_Mod/Form_Pressure_Matrix.f90"
#   include "Process_Mod/Insert_Volume_Source_For_Pressure.f90"

    ! Energy / Enthalpy
#   include "Process_Mod/Compute_Energy.f90"
#   include "Process_Mod/Insert_Energy_Bc.f90"

  end module
