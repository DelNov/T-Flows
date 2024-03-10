!==============================================================================!
  real function Monin_Obukov_Thermal(Turb,                   &
                                     u_tan, wall_dist, z_o,  &
                                     t_p,  t_wall, grav)
!------------------------------------------------------------------------------!
!   Calculate wall-heat flux according to Monin-Obukov Similarity              !
!   Theory (MOST)                                                              !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Turb_Type) :: Turb
  real, intent(in) :: u_tan
  real, intent(in) :: wall_dist
  real, intent(in) :: z_o
  real, intent(in) :: t_p
  real, intent(in) :: t_wall
  real, intent(in) :: grav
!------------------------------[Local parameters]------------------------------!
  real, parameter :: B1      = 9.4    ! parameters from Uno1995
  real, parameter :: B2      = 4.7
  real, parameter :: DTH     = 5.3
  real, parameter :: PR_TURB = 0.74   ! Businger et al. (1971)
!-----------------------------------[Locals]-----------------------------------!
  real :: ri_bo, fth, cth
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Turb)
!==============================================================================!

  ri_bo = grav*wall_dist*(t_p - t_wall) &
              / ((t_wall + 273.15) * u_tan**2)

  if((t_p - t_wall) > 0.0) then           ! stable condition
    fth = 1./(1. + B2*ri_bo)**2
  else if((t_p - t_wall) < 0.0) then      ! unstable condition
    cth = (DTH*kappa**2)/((log(wall_dist/z_o))**2)*B1*sqrt(wall_dist/z_o)
    fth = 1. - (DTH*ri_bo)/(1. + cth*sqrt(abs(ri_bo)))
  else if(abs(t_p - t_wall) < TINY) then  ! neutral condition
    fth = 1.0
  end if

  ! Specify the return value
  Monin_Obukov_Thermal = fth / PR_TURB

  end function
