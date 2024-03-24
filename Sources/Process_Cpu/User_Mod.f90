#include "../Shared/Assert.h90"
#include "../Shared/Browse.h90"

!==============================================================================!
  module User_Mod
!------------------------------------------------------------------------------!
!   This is embrio of a future User module, a place where user can
!   define his/her variables and pass them around his functions
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Porosity_Mod
  use Swarm_Mod
  use Interface_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

# include "User_Mod/Types.f90"

  contains

#   include "User_Mod/Before_Exit.f90"
#   include "User_Mod/Beginning_Of_Compute_Energy.f90"
#   include "User_Mod/Beginning_Of_Compute_Momentum.f90"
#   include "User_Mod/Beginning_Of_Compute_Pressure.f90"
#   include "User_Mod/Beginning_Of_Compute_Scalar.f90"
#   include "User_Mod/Beginning_Of_Correct_Velocity.f90"
#   include "User_Mod/Beginning_Of_Iteration.f90"
#   include "User_Mod/Beginning_Of_Simulation.f90"
#   include "User_Mod/Beginning_Of_Time_Step.f90"
#   include "User_Mod/Bulk_Velocity.f90"
#   include "User_Mod/Calculate_Mean.f90"
#   include "User_Mod/End_Of_Compute_Energy.f90"
#   include "User_Mod/End_Of_Compute_Momentum.f90"
#   include "User_Mod/End_Of_Compute_Pressure.f90"
#   include "User_Mod/End_Of_Compute_Scalar.f90"
#   include "User_Mod/End_Of_Correct_Velocity.f90"
#   include "User_Mod/End_Of_Iteration.f90"
#   include "User_Mod/End_Of_Simulation.f90"
#   include "User_Mod/End_Of_Time_Step.f90"
#   include "User_Mod/Initialize_Variables.f90"
#   include "User_Mod/Insert_Particles.f90"
#   include "User_Mod/Interface_Exchange.f90"
#   include "User_Mod/Force.f90"
#   include "User_Mod/Save_Results.f90"
#   include "User_Mod/Save_Swarm.f90"
#   include "User_Mod/Source.f90"

  end module


