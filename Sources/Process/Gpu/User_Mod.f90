#include "../../Shared/Assert.h90"
#include "../../Shared/Browse.h90"
#include "../../Shared/Macros.h90"

!==============================================================================!
  module User_Mod
!------------------------------------------------------------------------------!
!   This is embrio of a future User module, a place where user can
!   define his/her variables and pass them around his functions
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Time_Mod
  use Turb_Mod
  use Gpu_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

# include "User_Mod/Types.f90"

  contains

#   include "User_Mod/Initialize_Variables.f90"
#   include "User_Mod/End_Of_Simulation.f90"
#   include "User_Mod/Save_Results.f90"
#   include "User_Mod/Source.f90"

  end module


