!==============================================================================!
  real function Roughness_Coefficient(Turb, z_o_function)
!------------------------------------------------------------------------------!
!   Set lower limit to roughness coefficient based on wall distance.           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Turb_Type) :: Turb
  real             :: z_o_function
!==============================================================================!

  Roughness_Coefficient = Turb % z_o

  if(z_o_function > -TINY) then
    Roughness_Coefficient = z_o_function
  end if

  end function
