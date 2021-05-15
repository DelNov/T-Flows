!==============================================================================!
  subroutine Get_Gas_And_Liquid_Phase(Vof, g, l)
!------------------------------------------------------------------------------!
!   Names indices for gas and liquid to input arguments                        !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Vof_Type), target :: Vof
  integer                 :: g, l
!==============================================================================!

  ! Distinguish between liquid and vapor
  l = 1; g = 2
  if(Vof % phase_dens(g) > Vof % phase_dens(l)) then
    l = 2; g = 1
  end if

  end subroutine
