!==============================================================================!
  real function Tau_Wall_Rough_Walls(Turb, dens, u_tau, u_tan, wall_dist, z_o)
!------------------------------------------------------------------------------!
!   Calculates tau_wall for rough walls.                                       !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Turb_Type) :: Turb
  real             :: dens, u_tau, u_tan, wall_dist, z_o
!==============================================================================!

  Tau_Wall_Rough_Walls = dens * kappa * u_tau * u_tan  &
                       / log(((wall_dist + z_o) / z_o))

  end function
