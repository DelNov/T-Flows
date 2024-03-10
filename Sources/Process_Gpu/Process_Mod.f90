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
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !------------------!
  !   Process type   !
  !------------------!
  type Process_Type

    contains

      ! Related to momentum conservation
      procedure :: Compute_Momentum
      procedure :: Form_Diffusion_Matrix
      procedure :: Insert_Diffusion_Bc
      procedure :: Add_Advection_Term
      procedure :: Add_Inertial_Term
      procedure :: Add_Pressure_Term
      procedure :: Correct_Velocity

      ! Related to pressure solution
      procedure :: Compute_Pressure
      procedure :: Form_Pressure_Matrix
      procedure :: Insert_Volume_Source_For_Pressure

  end type

  type(Process_Type) :: Process

  contains

    ! Related to momentum conservation
#   include "Process_Mod/Compute_Momentum.f90"
#   include "Process_Mod/Form_Diffusion_Matrix.f90"
#   include "Process_Mod/Insert_Diffusion_Bc.f90"
#   include "Process_Mod/Add_Advection_Term.f90"
#   include "Process_Mod/Add_Inertial_Term.f90"
#   include "Process_Mod/Add_Pressure_Term.f90"
#   include "Process_Mod/Correct_Velocity.f90"

    ! Related to pressure solution
#   include "Process_Mod/Compute_Pressure.f90"
#   include "Process_Mod/Form_Pressure_Matrix.f90"
#   include "Process_Mod/Insert_Volume_Source_For_Pressure.f90"

  end module
