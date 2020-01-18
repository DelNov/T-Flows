!==============================================================================!
  module Save_Results_Mod
!------------------------------------------------------------------------------!
!   Module containig functions for saving numerical results for visualization. !
!   It comes in two flavors: "Vtu" and "Cgns", depending on the file format    !
!   one wants to save.  It has (and uses) a sister module "Save_Grid_Mod",     !
!   which is in the directory "Shared".                                        !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Comm_Mod
  use Cpu_Timer_Mod,  only: Cpu_Timer_Mod_Start, Cpu_Timer_Mod_Stop
  use File_Mod
  use Grid_Mod,       only: Grid_Type
  use Cgns_Mod
  use Save_Grid_Mod
  use Turb_Mod,       NO_TURBULENCE => NONE
  use Swarm_Mod,      only: Particle_Type, Swarm_Type
  use Surf_Mod,       only: Vert_Type, Elem_Type, Surf_Type
  use Multiphase_Mod, NO_MULTIPHASE => NONE
  use User_Mod, only: n_user_arrays, user_array
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  contains

  include 'Save_Results_Mod/Cgns/Save_Results.f90'
  include 'Save_Results_Mod/Cgns/Save_Swarm.f90'
  include 'Save_Results_Mod/Cgns/Save_Surf.f90'

  end module 
