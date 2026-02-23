!==============================================================================!
  pure real function U_Tau_Log_Law(Turb, u_tan, wall_dist, nu, z_o)
!------------------------------------------------------------------------------!
!   Calculates u_tau from log law using Newton iterations (smooth wall)        !
!   and explicit formula for rough wall with z_o.                              !
!
!   Smooth wall equation:
!     u_tan / u_tau = (1/kappa) * ln( E * wall_dist * u_tau / nu )
!
!   Rough wall equation used consistently with your U_Plus_Log_Law:
!     U+ = (1/kappa) ln( (wall_dist + z_o) / z_o )
!     => u_tau = u_tan / U+
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Turb_Type), intent(in) :: Turb
  real,             intent(in) :: u_tan, wall_dist, nu, z_o
!---------------------------------[Locals]-------------------------------------!
  real    :: utau, utau_new, invk, arg, f, denom, uplus
  integer :: it
  real, parameter :: REL_TOL = 1.0e-6
  integer, parameter :: IT_MAX = 20
!==============================================================================!

  invk = 1.0 / (Turb % kappa + TINY)

  ! Trivial case
  if(u_tan <= TINY .or. wall_dist <= TINY .or. nu <= TINY) then
    U_Tau_Log_Law = 0.0
    return
  end if

  ! Rough wall: explicit
  if(z_o > TINY) then
    uplus = log( (wall_dist + z_o) / z_o ) * invk + TINY
    U_Tau_Log_Law = max(u_tan / max(uplus, TINY), 0.0)
    return
  end if

  ! Smooth wall: Newton iterations
  ! Initial guess (možeš zamijeniti svojom power-law procjenom ako želiš)
  utau = max(1.0e-8, 0.05 * u_tan)

  do it = 1, IT_MAX

    ! arg = E * y+  = E * wall_dist * utau / nu
    arg = Turb % e_log * wall_dist * utau / nu

    ! mimic your U_Plus_Log_Law safeguard: y_plus >= 1.05
    arg = max(arg, 1.05 * Turb % e_log)

    ! F(utau) = u_tan/utau - (1/kappa) ln(arg)
    f = u_tan/utau - invk * log(arg)

    ! stable Newton denominator: (u_tan/utau) + (1/kappa)
    denom = (u_tan/utau) + invk

    utau_new = utau + utau * f / max(denom, TINY)

    ! Safety
    utau_new = max(utau_new, 1.0e-12)

    if(abs(utau_new - utau) / max(utau, 1.0e-12) < REL_TOL) then
      utau = utau_new
      exit
    end if

    utau = utau_new
  end do

  U_Tau_Log_Law = utau

  end function

