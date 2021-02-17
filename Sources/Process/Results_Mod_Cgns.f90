!==============================================================================!
  module Results_Mod
!------------------------------------------------------------------------------!
!   Module containig functions for saving numerical results for visualization. !
!   It comes in two flavors: "Vtu" and "Cgns", depending on the file format    !
!   one wants to save.  It has (and uses) a sister module "Save_Grid_Mod",     !
!   which is in the directory "Shared".                                        !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Save_Grid_Mod
  use Backup_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !------------------!
  !   Results type   !
  !------------------!
  type Results_Type
    integer :: interval  ! result save interval
    logical :: initial   ! save intial condition or not
  end type
  type(Results_Type) :: result

  contains

  include 'Results_Mod/Main.f90'
  include 'Results_Mod/Time_To_Save.f90'
  include 'Results_Mod/Cgns/Save.f90'
  include 'Results_Mod/Cgns/Save_Swarm.f90'
  include 'Results_Mod/Cgns/Save_Surf.f90'

  end module 
