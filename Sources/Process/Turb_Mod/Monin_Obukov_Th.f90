!==============================================================================!
  real function Monin_Obukov_Th(Turb, c1, c2, u_tan, wall_dist, z_o, t_p, & 
                  t_wall, grav)
!------------------------------------------------------------------------------!
!   Calculate wall-heat flux according to Monin-Obukov Similarity              !
!   Theory (MOST)                                                              !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Turb_Type)    :: Turb
  integer, intent(in) :: c1, c2
!-----------------------------------[Locals]-----------------------------------!
  real :: z_o, Ri_bo, Fth, cth, wall_dist, t_p, t_wall, u_tan, grav
  real, parameter :: b1 = 9.4          ! parameters from Uno1995
  real, parameter :: b2 = 4.7
  real, parameter :: dth = 5.3
  real, parameter :: Pr_turb = 0.74    ! Businger et al. (1971)
  real, parameter :: prandtlmol = 0.71
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Turb)
!==============================================================================!

  ! Take the value specified in control file

  Ri_bo = grav*wall_dist*(t_p - t_wall) &
              / ((t_wall + 273.15) * u_tan**2)

  if((t_p - t_wall) > 0.0) then           ! Stable condition
    Fth = 1./(1. + b2*Ri_bo)**2
  else if((t_p - t_wall) < 0.0) then      ! Unstable condition
    cth = (dth*kappa**2)/((log(wall_dist/z_o))**2)*b1*sqrt(wall_dist/z_o)
    Fth = 1. - (dth*Ri_bo)/(1. + cth*sqrt(abs(Ri_bo)))
  else if(abs(t_p - t_wall) < tiny) then  ! Neutral condition
    Fth = 1.0
  end if

  ! Specify the return value
  Monin_Obukov_Th = (kappa**2/(LOG(wall_dist/z_o))**2 * (t_p - t_wall) &
                   * u_tan * Fth/Pr_turb)

  end function
