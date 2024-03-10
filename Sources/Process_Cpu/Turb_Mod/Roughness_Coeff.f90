!==============================================================================!
  real function Roughness_Coeff(Turb, c1, c2)
!------------------------------------------------------------------------------!
!   Set up roughness coefficient                                               !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Turb_Type)    :: Turb
  integer, intent(in) :: c1, c2
!-----------------------------------[Locals]-----------------------------------!
  real :: z_o
!==============================================================================!

  ! Take the value specified in control file
  z_o = Turb % z_o(c2)

  ! Set lower limit to roughness coefficient based on wall distance
  if(z_o .gt. TINY) then
    z_o = max(Turb % pnt_grid % wall_dist(c1)  &
        / (e_log * max(Turb % y_plus(c1), 1.0)), z_o)
  end if

  ! Specify the return value
  Roughness_Coeff = z_o

  end function
