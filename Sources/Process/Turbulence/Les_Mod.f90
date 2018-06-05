!==============================================================================!
  module Les_Mod
!------------------------------------------------------------------------------!
!   Definition of variables used by LES turbulence models.                     !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Var_Mod
  use Turbulence_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  ! Variables relevant for 'LES' computations
  real              :: c_smag
  real, allocatable :: c_dyn(:)

  ! For LES you need to know nearest wall cell
  integer, allocatable :: nearest_wall_cell(:)

! real, allocatable :: kin_sgs(:)
! real, allocatable :: vis_t_sgs(:)

  real, allocatable :: wale_v(:)

end module  
