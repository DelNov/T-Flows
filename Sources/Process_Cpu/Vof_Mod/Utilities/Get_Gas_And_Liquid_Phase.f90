!==============================================================================!
  subroutine Get_Vapour_And_Liquid_Phase(Vof, vapour, liquid)
!------------------------------------------------------------------------------!
!   Distinguish between liquid and vapor phases, based on their density        !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Vof_Type), target :: Vof
  integer                 :: vapour, liquid
!==============================================================================!

  ! Assume vapour is zero and liquid is one ...
  vapour = 0; liquid = 1

  ! ... and if you happened to make a wrong assumption, swap them around
  if(Vof % phase_dens(vapour) > Vof % phase_dens(liquid)) then
    call Swap_Int(vapour, liquid)
  end if

  end subroutine
