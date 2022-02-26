!==============================================================================!
  real function Tau_Wall_Low_Re(Turb, dens, u_tau, u_tan, y_plus)
!------------------------------------------------------------------------------!
!   Calculates tau_wall for low Reynolds approach.                             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Turb_Type) :: Turb
  real             :: dens, u_tau, u_tan, y_plus
!==============================================================================!

  Tau_Wall_Low_Re = dens * kappa * u_tau * u_tan   &
                  / log(e_log * max(y_plus, 1.05))

  end function
