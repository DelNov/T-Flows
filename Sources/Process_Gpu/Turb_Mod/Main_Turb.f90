!==============================================================================!
  subroutine Main_Turb(Turb, Flow, Grid)
!------------------------------------------------------------------------------!
!   Turbulence model main function (called inside inner iterations)            !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Turb_Type) :: Turb  !! parent class
  type(Field_Type) :: Flow  !! flow field object
  type(Grid_Type)  :: Grid  !! grid object
!==============================================================================!

  !---------------------------------------------------!
  !   Start branching for various turbulence models   !
  !---------------------------------------------------!

  if(Turb % model .eq. LES_SMAGORINSKY) then
    call Flow % Calculate_Shear_And_Vorticity(Grid)
  end if

  end subroutine
