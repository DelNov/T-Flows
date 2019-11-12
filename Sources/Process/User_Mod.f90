!==============================================================================!
  module User_Mod
!------------------------------------------------------------------------------!
!   This is embrio of a future User module, a place where user can
!   define his/her variables and pass them around his functions
!------------------------------------------------------------------------------!
  use Const_Mod
  use File_Mod
  use Comm_Mod
  use Grid_Mod,    only: Grid_Type
  use Matrix_Mod,  only: Matrix_Type
  use Var_Mod,     only: Var_Type
  use Field_Mod
  use Bulk_Mod,    only: Bulk_Type
  use Turb_Mod                       ! needed for saving results
  use Swarm_Mod
  use Multiphase_Mod
  use Control_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  integer           :: n_user_arrays
  real, allocatable :: user_array(:,:)

  contains

  include 'User_Mod/Allocate.f90'
  include 'User_Mod/Before_Exit.f90'
  include 'User_Mod/Beginning_Of_Time_Step.f90'
  include 'User_Mod/Beginning_Of_Compute_Energy.f90'
  include 'User_Mod/Beginning_Of_Compute_Momentum.f90'
  include 'User_Mod/Beginning_Of_Compute_Pressure.f90'
  include 'User_Mod/Beginning_Of_Correct_Velocity.f90'
  include 'User_Mod/Calculate_Mean.f90'
  include 'User_Mod/End_Of_Compute_Energy.f90'
  include 'User_Mod/End_Of_Compute_Momentum.f90'
  include 'User_Mod/End_Of_Compute_Pressure.f90'
  include 'User_Mod/End_Of_Correct_Velocity.f90'
  include 'User_Mod/End_Of_Time_Step.f90'
  include 'User_Mod/Initialize.f90'
  include 'User_Mod/Force.f90'
  include 'User_Mod/Save_Results.f90'
  include 'User_Mod/Source.f90'

  end module 


