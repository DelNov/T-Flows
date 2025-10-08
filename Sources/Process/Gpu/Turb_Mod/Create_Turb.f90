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

  ! Create deltas
  if(Turb % model .eq. SPALART_ALLMARAS .or.  &
     Turb % model .eq. DES_SPALART) then
    allocate(Turb % h_max(-nb:nc));  Turb % h_max = 0.
    allocate(Turb % h_min(-nb:nc));  Turb % h_min = 0.
    allocate(Turb % h_w  (-nb:nc));  Turb % h_w   = 0.
  end if

  ! Hydraulic roughness
  allocate(Turb % z_o(-nb:nc)); Turb % z_o = 0.0

  !----------------------!
  !   Spalart Allmaras   !
  !----------------------!
  if(Turb % model .eq. SPALART_ALLMARAS .or.  &
     Turb % model .eq. DES_SPALART) then

    call Var_Mod_Create_Variable(Turb % vis, Grid, 'VIS', '')

    ! Other variables such as time scale, length scale and production
    allocate(Turb % vis_t(-nb:nc));  Turb % vis_t = 0.
    allocate(Turb % vis_w(-nb:nc));  Turb % vis_w = 0.  ! wall visc

    if(Flow % heat_transfer) then
      allocate(Turb % con_w(-nb:nc));  Turb % con_w = 0.  ! wall cond
    end if ! Flow % heat_transfer

  end if

  !-----------------------!
  !   Smagorinsky model   !
  !-----------------------!
  if(Turb % model .eq. LES_SMAGORINSKY) then

    ! Variables such as time scale, length scale and production
    allocate(Turb % vis_t(-nb:nc));  Turb % vis_t = 0.
    allocate(Turb % vis_w(-nb:nc));  Turb % vis_w = 0.

    if(Flow % heat_transfer) then
      allocate(Turb % con_w(-nb:nc));  Turb % con_w = 0.
    end if

  end if ! LES_SMAGORINSKY

  !----------------!
  !   Wale model   !
  !----------------!
  if(Turb % model .eq. LES_WALE) then

    allocate(Turb % wale_v(-nb:nc));  Turb % wale_v = 0.

    ! Other variables such as time scale, length scale and production
    allocate(Turb % vis_t   (-nb:nc));  Turb % vis_t = 0.
    allocate(Turb % vis_w   (-nb:nc));  Turb % vis_w = 0.

    if(Flow % heat_transfer) then
      allocate(Turb % con_w(-nb:nc));  Turb % con_w = 0.
    end if

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
