!==============================================================================!
  real function Tau_Wall_Log_Law(Turb, dens, u_tau, u_tan, wall_dist,  &
                                 y_plus, z_o)
!------------------------------------------------------------------------------!
!   Calculates wall shear stress by using log law.                             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Turb_Type) :: Turb
  real             :: dens, u_tau, u_tan, wall_dist, y_plus, z_o
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Turb)
!==============================================================================!

  if(z_o > TINY) then

    Tau_Wall_Log_Law = dens * kappa * u_tau * u_tan  &
                     / log(((wall_dist + z_o) / z_o))
  else

    Tau_Wall_Log_Law = dens * kappa * u_tau * u_tan   &
                     / log(e_log * max(y_plus, 1.05))
  end if

  end function
