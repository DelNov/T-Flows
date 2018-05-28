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
  real, allocatable :: c_dyn(:), c_dyn_mean(:)

  ! For LES you need to know nearest wall cell
  integer, allocatable :: nearest_wall_cell(:)

  real, allocatable :: shear_mean(:), kin_sgs(:)
  real, allocatable :: vis_t_sgs(:), vis_t_mean(:)

  real, allocatable :: shear_r(:), shear_mean_r(:), wale_v(:)

end module  
