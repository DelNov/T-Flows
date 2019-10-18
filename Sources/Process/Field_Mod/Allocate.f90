!==============================================================================!
  subroutine Field_Mod_Allocate(flow, grid)
!------------------------------------------------------------------------------!
!   Allocates memory for the entire field.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type)         :: flow
  type(Grid_Type),  target :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer          :: sc
  character(len=4) :: c_name, q_name
!==============================================================================!

  ! Store the pointer to a grid
  flow % pnt_grid => grid

  !----------------------------!
  !   Navier-Stokes equation   !
  !----------------------------!

  ! Allocate memory for velocity components
  call Var_Mod_Allocate_Solution('U', '', flow % u, grid)
  call Var_Mod_Allocate_Solution('V', '', flow % v, grid)
  call Var_Mod_Allocate_Solution('W', '', flow % w, grid)

  ! Allocate memory for pressure correction and pressure
  call Var_Mod_Allocate_New_Only('PP', flow % pp, grid)
  call Var_Mod_Allocate_New_Only('P',  flow % p,  grid)

  ! Allocate memory for volumetric flux
  call Var_Mod_Allocate_Vol_Flux('V_FLUX','', flow % vol_flux, grid)

  ! Mass flow rates at cell faces are always needed
  allocate(flow % flux(grid % n_faces));  flow % flux = 0.

  ! density at cell faces:
  allocate(dens_face(grid % n_faces))

  !-----------------------------------------!
  !   Enthalpy conservation (temperature)   !
  !-----------------------------------------!
  if(heat_transfer) then
    call Var_Mod_Allocate_Solution('T', 'Q', flow % t, grid)
  end if ! heat_transfer

  allocate(flow % vort (-grid % n_bnd_cells:grid % n_cells)); flow % vort  = 0.
  allocate(flow % shear(-grid % n_bnd_cells:grid % n_cells)); flow % shear = 0.

  !--------------------------------------!
  !   Allocate memory for user scalars   !
  !--------------------------------------!
  allocate(flow % scalar(flow % n_scalars))

  !-------------------------------------!
  !   Browse through all user scalars   !
  !-------------------------------------!
  do sc = 1, flow % n_scalars

    ! Set variable name
    c_name = 'C_00'
    q_name = 'Q_00'
    write(c_name(3:4),'(i2.2)') sc
    write(q_name(3:4),'(i2.2)') sc

    ! Allocate memory for passive scalar
    call Var_Mod_Allocate_Solution(c_name, q_name, flow % scalar(sc), grid)

  end do

  end subroutine
