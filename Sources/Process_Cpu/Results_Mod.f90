#include "../Shared/Assert.h90"
#include "../Shared/Browse.h90"
#include "../Shared/Unused.h90"

!==============================================================================!
  module Results_Mod
!------------------------------------------------------------------------------!
!>  The Results_Mod module in T-Flows is dedicated to managing and saving
!>  numerical results for visualization and backup (needed for restart).
!>  It defines the Results_Type which encapsulates data members used to define
!>  the frequency of saving and wheather a separate set of results will be
!>  saved for boundary and a number of member functions to control the logic
!>  of saving or perform some particular saving tasks.  The module was creted
!>  to minimize the complexity within the main function of the T-Flows'
!>  sub-program Process.
!------------------------------------------------------------------------------!
!   Key Features                                                               !
!                                                                              !
!   * Control of result saving: Allows the user to specify parameters          !
!     controlling the frequency and conditions under which results are saved.  !
!     This includes options to save results at boundaries, during initial      !
!     conditions, at regular intervals, and at specific intervals for particle !
!     swarms.                                                                  !
!   * Member functions: A collection of procedures to handle different aspects !
!     of result saving. These functions are responsible for saving results in  !
!     .vtu format, handling scalar, vector, and tensor data, as well as        !
!     managing the output of front tracking and surface data.                  !
!   * Utilizes Backup_Mod for backing up and restoring simulation states.      !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Backup_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !------------------!
  !   Results type   !
  !------------------!
  !> Encapsulates data members to control the frequency and type of saving
  !> and procedures to control the saving process and different aspects of
  !> saving results.
  type Results_Type

    logical :: boundary        !! set to true to save results at boundaries
    logical :: units           !! set to true to save variable name with unit
    logical :: initial         !! set to true to save intial condition
    integer :: interval        !! result save interval
    integer :: interval_swarm  !! result save interval for particles

    contains
      procedure :: Main_Results

      procedure, private :: Save_Vtu_Fields
      procedure, private :: Save_Vtu_Front
      procedure, private :: Save_Vtu_Header_Int
      procedure, private :: Save_Vtu_Scalar_Int
      procedure, private :: Save_Vtu_Scalar_Real
      procedure, private :: Save_Vtu_Surf
      procedure, private :: Save_Vtu_Swarm
      procedure, private :: Save_Vtu_Tensor_6_Real
      procedure, private :: Save_Vtu_Tensor_9_Real
      procedure, private :: Save_Vtu_Vector_Real
      procedure, private :: Time_To_Save_Results
      procedure, private :: Time_To_Save_Swarm
      procedure, private :: Var_Name

  end type

  type(Results_Type) :: Results

  contains

#   include "Results_Mod/Main_Results.f90"
#   include "Results_Mod/Save_Vtu_Fields.f90"
#   include "Results_Mod/Save_Vtu_Front.f90"
#   include "Results_Mod/Save_Vtu_Header_Int.f90"
#   include "Results_Mod/Save_Vtu_Scalar_Int.f90"
#   include "Results_Mod/Save_Vtu_Scalar_Real.f90"
#   include "Results_Mod/Save_Vtu_Surf.f90"
#   include "Results_Mod/Save_Vtu_Swarm.f90"
#   include "Results_Mod/Save_Vtu_Tensor_6_Real.f90"
#   include "Results_Mod/Save_Vtu_Tensor_9_Real.f90"
#   include "Results_Mod/Save_Vtu_Vector_Real.f90"
#   include "Results_Mod/Time_To_Save_Results.f90"
#   include "Results_Mod/Time_To_Save_Swarm.f90"
#   include "Results_Mod/Var_Name.f90"

  end module
