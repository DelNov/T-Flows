!==============================================================================!
  module Results_Mod
!------------------------------------------------------------------------------!
!   Module containig functions for saving numerical results for visualization. !
!   It comes in two flavors: "Vtu" and "Cgns", depending on the file format    !
!   one wants to save.  It has (and uses) a sister module "Save_Grid_Mod",     !
!   which is in the directory "Shared".                                        !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Cpu_Timer_Mod,  only: Cpu_Timer_Mod_Start, Cpu_Timer_Mod_Stop
  use File_Mod
  use Save_Grid_Mod
  use Swarm_Mod,      only: Particle_Type, Swarm_Type
  use Multiphase_Mod
  use Backup_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !------------------!
  !   Results type   !
  !------------------!
  type Results_Type
    integer :: interval
  end type
  type(Results_Type) :: result

  contains

  include 'Results_Mod/Time_To_Save.f90'
  include 'Results_Mod/Vtu/Save.f90'              ! binary
  include 'Results_Mod/Vtu/Save_Scalar_Int.f90'   ! binary
  include 'Results_Mod/Vtu/Save_Scalar_Real.f90'  ! binary
  include 'Results_Mod/Vtu/Save_Swarm.f90'
  include 'Results_Mod/Vtu/Save_Vector_Real.f90'  ! binary
  include 'Results_Mod/Vtu/Save_Surf.f90'

  end module 
