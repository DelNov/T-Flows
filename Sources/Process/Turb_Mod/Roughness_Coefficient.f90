!==============================================================================!
  real function Roughness_Coefficient(turb, z_o_function)
!------------------------------------------------------------------------------!
!   Set lower limit to roughness coefficient based on wall distance.           !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Turb_Mod,  only: Turb_Type
  use Const_Mod, only: TINY
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Turb_Type) :: turb
  real            :: z_o_function
!==============================================================================!

  Roughness_Coefficient = turb % z_o

  if(z_o_function > -TINY) then
    Roughness_Coefficient = z_o_function
  end if

  end function
