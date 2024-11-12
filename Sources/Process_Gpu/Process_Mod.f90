#include "../Shared/Assert.h90"
#include "../Shared/Browse.h90"
#include "../Shared/Macros.h90"
#include "../Shared/Unused.h90"

!==============================================================================!
  module Process_Mod
!------------------------------------------------------------------------------!
  use Assert_Mod
  use Profiler_Mod
  use Iter_Mod
  use Gpu_Mod
  use Control_Mod
  use Turb_Mod
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

      ! Splash screen
      procedure :: Logo_Pro

      ! Related to momentum conservation
      procedure :: Compute_Momentum
      procedure :: Form_Momentum_Matrix
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
      procedure :: Form_Energy_Matrix
      procedure :: Insert_Energy_Bc

  end type

  type(Process_Type) :: Process  !! singleton object of Process_Type, introduced
                                 !! to allow easy access to its procedures
  contains

    ! Splash screen
#   include "Process_Mod/Logo_Pro.f90"

    ! Related to momentum conservation
#   include "Process_Mod/Compute_Momentum.f90"
#   include "Process_Mod/Form_Momentum_Matrix.f90"
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
#   include "Process_Mod/Form_Energy_Matrix.f90"
#   include "Process_Mod/Insert_Energy_Bc.f90"

  end module
