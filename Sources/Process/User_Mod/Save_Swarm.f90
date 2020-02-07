!==============================================================================!
   subroutine User_Mod_Save_Swarm(flow, turb, mult, swarm, n) 
!------------------------------------------------------------------------------!
!    User's Save_Swarm function.                                               !
!------------------------------------------------------------------------------!
  use Const_Mod                      ! constants
  use Comm_Mod                       ! parallel stuff
  use Grid_Mod,  only: Grid_Type
  use Field_Mod, only: Field_Type
  use Bulk_Mod,  only: Bulk_Type
  use Var_Mod,   only: Var_Type
  use File_Mod,  only: problem_name
  use Turb_Mod 
  use Swarm_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),       target  :: flow
  type(Turb_Type),        target  :: turb
  type(Multiphase_Type),  target  :: mult
  type(Swarm_Type),       target  :: swarm
  integer                         :: n      ! current time step
!-----------------------------------[Locals]-----------------------------------!
  type(Var_Type),  pointer  :: u, v, w, t
  type(Grid_Type), pointer  :: grid
  type(Bulk_Type), pointer  :: bulk
!==============================================================================!

  ! Take aliases for the flow 
  grid => flow % pnt_grid
  bulk => flow % bulk
  call Field_Mod_Alias_Momentum(flow, u, v, w)
  call Field_Mod_Alias_Energy  (flow, t)

  end subroutine
