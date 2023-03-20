!==============================================================================!
  real function Y_Plus_Rough_Walls(Turb, u_tau, wall_dist, kin_vis, z_o)
!------------------------------------------------------------------------------!
!   Calculates y+ for rough walls.                                             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Turb_Type) :: Turb
  real             :: u_tau, wall_dist, kin_vis, z_o
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Turb)
!==============================================================================!

  Y_Plus_Rough_Walls = u_tau * (wall_dist + z_o) / kin_vis

  end function
