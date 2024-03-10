!==============================================================================!
  real function U_Plus_Log_Law(Turb, wall_dist, y_plus, z_o)
!------------------------------------------------------------------------------!
!   Calculates U+ from log law.                                                !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Turb_Type) :: Turb
  real             :: y_plus, wall_dist, z_o
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Turb)
!==============================================================================!

  if(z_o > TINY) then

    U_Plus_Log_Law = log( (wall_dist + z_o) / z_o)  &
                      / (kappa + TINY) + TINY
  else

    U_Plus_Log_Law = log( max(y_plus, 1.05) * e_log ) / kappa
  end if

  end function
