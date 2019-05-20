!==============================================================================!
  real function Roughness_Coefficient(z_o_constant, z_o_function, c)
!------------------------------------------------------------------------------!
!   Set lower limit to roughness coefficient based on wall distance.           !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod, only: TINY
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real    :: z_o_constant, z_o_function
  integer :: c
!==============================================================================!

  Roughness_Coefficient = z_o_constant

  if(z_o_function > TINY) then
    Roughness_Coefficient = z_o_function
  end if

  end function
