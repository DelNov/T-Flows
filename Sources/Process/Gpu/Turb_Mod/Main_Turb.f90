!==============================================================================!
  subroutine Main_Turb(Turb, Grid, Flow)
!------------------------------------------------------------------------------!
!   Turbulence model main function (called inside inner iterations)            !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Turb_Type) :: Turb  !! parent class
  type(Grid_Type)  :: Grid  !! grid object
  type(Field_Type) :: Flow  !! flow field object
!==============================================================================!

  !---------------------------------------------------!
  !   Start branching for various turbulence models   !
  !---------------------------------------------------!

  if(Turb % model .eq. LES_SMAGORINSKY .or.  &
     Turb % model .eq. LES_WALE) then
    call Flow % Calculate_Shear_And_Vorticity(Grid)
  end if

  if(Turb % model .eq. SPALART_ALLMARAS .or.  &
     Turb % model .eq. DES_SPALART) then
    call Flow % Calculate_Shear_And_Vorticity(Grid)

    call Turb % Compute_Variable(Grid, Turb % vis, Flow)
    call Turb % Vis_T_Spalart_Allmaras(Grid, Flow)
  end if

  end subroutine
