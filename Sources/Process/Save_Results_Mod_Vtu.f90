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
  use Save_Grid_Mod
  use Turb_Mod
  use Swarm_Mod,      only: Particle_Type, Swarm_Type
  use Surf_Mod,       only: Vert_Type, Elem_Type, Surf_Type
  use Multiphase_Mod
  use User_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  contains

  include 'Save_Results_Mod/Vtu/Save_Results.f90'      ! binary
  include 'Save_Results_Mod/Vtu/Save_Scalar_Int.f90'   ! binary
  include 'Save_Results_Mod/Vtu/Save_Scalar_Real.f90'  ! binary
  include 'Save_Results_Mod/Vtu/Save_Swarm.f90'
  include 'Save_Results_Mod/Vtu/Save_Vector_Real.f90'  ! binary
  include 'Save_Results_Mod/Vtu/Save_Surf.f90'

  end module 
