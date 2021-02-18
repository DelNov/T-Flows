!==============================================================================!
  subroutine Multiphase_Mod_Vof_Get_Gas_And_Liquid_Phase(mult, g, l)
!------------------------------------------------------------------------------!
!   Names indices for gas and liquid to input arguments                        !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Multiphase_Type), target :: mult
  integer                       :: g, l
!==============================================================================!

  ! Distinguish between liquid and vapor
  l = 1; g = 2
  if(mult % phase_dens(g) > mult % phase_dens(l)) then
    l = 2; g = 1
  end if

  end subroutine
