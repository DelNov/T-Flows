!==============================================================================!
  subroutine Init_Turb(Turb, Grid, Flow)
!------------------------------------------------------------------------------!
!>  Turbulence model initializations (at the beginning of a time step)         !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Turb_Type) :: Turb  !! parent class
  type(Grid_Type)  :: Grid  !! numerical grid
  type(Field_Type) :: Flow  !! flow field
!==============================================================================!

  !--------------------------------------------------------!
  !   Start initializations of various turbulence models   !
  !--------------------------------------------------------!

  if(Turb % model .eq. LES_SMAGORINSKY .or.  &
     Turb % model .eq. LES_WALE) then
    call Flow % Calculate_Shear_And_Vorticity(Grid)
    if(Turb % model .eq. LES_WALE) then
      call Turb % Vis_T_Wale(Flow, Grid)
    end if
    call Turb % Vis_T_Subgrid(Flow, Grid)
  end if

  end subroutine
