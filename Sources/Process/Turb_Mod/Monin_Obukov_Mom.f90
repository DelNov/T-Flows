!==============================================================================!
  real function Monin_Obukov_Mom(Turb, c1, c2, u_tan, wall_dist, z_o, t_p, & 
                  t_wall, grav)
!------------------------------------------------------------------------------!
!   Calculate friction velocity u_tau according to Monin-Obukov Similarity     !
!   Theory (MOST)                                                              !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Turb_Type)    :: Turb
  integer, intent(in) :: c1, c2
!-----------------------------------[Locals]-----------------------------------!
  real :: z_o, Ri_bo, Fm, cm, wall_dist, t_p, t_wall, u_tan, grav
  real, parameter :: b1 = 9.4 !parameters from Uno1995
  real, parameter :: b2 = 4.7
  real, parameter :: dm = 7.4
  real, parameter :: prandtlmol = 0.71
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Turb)
!==============================================================================!

  Ri_bo = grav*wall_dist*(t_p - t_wall) &
          / ((t_wall + 273.15) * u_tan**2)

  if((t_p - t_wall) > 0.0) then               ! Stable condition
    Fm = 1./(1. + b2*Ri_bo)**2
  else if((t_p - t_wall) < 0.0) then          ! Unstable condition
    cm = (dm*kappa**2)/((log(wall_dist/z_o))**2)*b1*sqrt(wall_dist/z_o)
    Fm = 1. - (dm*Ri_bo)/(1. + cm*sqrt(abs(Ri_bo)))
  else if(abs(t_p - t_wall) < tiny) then      ! Neutral condition
    Fm = 1.0
  end if

  ! Specify the return value
  Monin_Obukov_Mom = (kappa/log(wall_dist/z_o) * u_tan * sqrt(Fm))**2

  end function
