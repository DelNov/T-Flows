!==============================================================================!
  subroutine Create_Turb(Turb, Grid, Flow)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Turb_Type), target :: Turb
  type(Grid_Type)          :: Grid
  type(Field_Type), target :: Flow
!-----------------------------------[Locals]-----------------------------------!
  integer :: nb, nc
!==============================================================================!

  ! Give some sign
  O_Print '(a)', ' # Creating the turbulence module'

  ! Take aliases
  nb = Grid % n_bnd_cells
  nc = Grid % n_cells

  !-----------------------!
  !   Smagorinsky model   !
  !-----------------------!
  if(Turb % model .eq. LES_SMAGORINSKY) then

    ! Variables such as time scale, length scale and production
    allocate(Turb % vis_t(-nb:nc));  Turb % vis_t = 0.

  end if ! LES_SMAGORINSKY

  !----------------!
  !   Wale model   !
  !----------------!
  if(Turb % model .eq. LES_WALE) then

    allocate(Turb % wale_v(-nb:nc));  Turb % wale_v = 0.

    ! Other variables such as time scale, length scale and production
    allocate(Turb % vis_t   (-nb:nc));  Turb % vis_t   = 0.

  end if

  !------------------------------------------------!
  !   Variables needed for all turbulence models   !
  !------------------------------------------------!
  if(Turb % model .ne. NO_TURBULENCE_MODEL) then

    allocate(Turb % y_plus(-nb:nc));  Turb % y_plus = 0.

    ! Shear and vorticity are defined in the Field_Mod
    allocate(Flow % shear(-nb:nc));  Flow % shear = 0.0
    allocate(Flow % vort (-nb:nc));  Flow % vort  = 0.0
  end if

  end subroutine
