#include "../../Shared/Assert.h90"
#include "../../Shared/Browse.h90"
#include "../../Shared/Unused.h90"

!==============================================================================!
  module Results_Mod
!----------------------------------[Modules]-----------------------------------!
  use Backup_Mod
  use User_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !------------------!
  !   Results type   !
  !------------------!
  type Results_Type

    logical :: boundary        !! set to true to save results at boundaries
    logical :: units           !! set to true to save variable name with unit
    logical :: initial         !! set to true to save intial condition
    integer :: interval        !! result save interval
    integer :: interval_swarm  !! result save interval for particles

    contains
      procedure :: Main_Results

      procedure, private :: Save_Vtu_Fields
      procedure, private :: Save_Vtu_Header_Int
      procedure, private :: Save_Vtu_Scalar_Int
      procedure, private :: Save_Vtu_Scalar_Real
      procedure, private :: Save_Vtu_Tensor_6_Real
      procedure, private :: Save_Vtu_Tensor_9_Real
      procedure, private :: Save_Vtu_Vector_Real
      procedure, private :: Time_To_Save_Results
      procedure, private :: Var_Name

  end type

  type(Results_Type) :: Results

  contains

#   include "Results_Mod/Main_Results.f90"
#   include "Results_Mod/Save_Vtu_Fields.f90"
#   include "Results_Mod/Save_Vtu_Header_Int.f90"
#   include "Results_Mod/Save_Vtu_Scalar_Int.f90"
#   include "Results_Mod/Save_Vtu_Scalar_Real.f90"
#   include "Results_Mod/Save_Vtu_Tensor_6_Real.f90"
#   include "Results_Mod/Save_Vtu_Tensor_9_Real.f90"
#   include "Results_Mod/Save_Vtu_Vector_Real.f90"
#   include "Results_Mod/Time_To_Save_Results.f90"
#   include "Results_Mod/Var_Name.f90"

  end module
