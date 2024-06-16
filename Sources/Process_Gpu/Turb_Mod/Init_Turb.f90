!==============================================================================!
  subroutine Init_Turb(Turb, Flow, Grid)
!------------------------------------------------------------------------------!
!>  Turbulence model initializations (at the beginning of a time step)         !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Turb_Type) :: Turb  !! parent class
  type(Field_Type) :: Flow  !! flow field
  type(Grid_Type)  :: Grid  !! numerical grid
!==============================================================================!

  !--------------------------------------------------------!
  !   Start initializations of various turbulence models   !
  !--------------------------------------------------------!

  if(Turb % model .eq. LES_SMAGORINSKY) then
    call Turb % Vis_T_Subgrid(Flow, Grid)
  end if

  end subroutine
