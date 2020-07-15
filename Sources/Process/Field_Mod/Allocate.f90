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

  call Field_Mod_Allocate_Grad_Matrix(flow)

  !----------------------------!
  !   Navier-Stokes equation   !
  !----------------------------!

  ! Allocate memory for velocity components
  call Var_Mod_Allocate_Solution(flow % u,   grid, 'U', '')
  call Var_Mod_Allocate_Solution(flow % v,   grid, 'V', '')
  call Var_Mod_Allocate_Solution(flow % w,   grid, 'W', '')

  ! Potential for initial velocity computation
  call Var_Mod_Allocate_Solution(flow % pot, grid, 'POT', '')

  ! Allocate memory for pressure correction and pressure
  call Var_Mod_Allocate_New_Only(flow % pp, grid, 'PP')
  call Var_Mod_Allocate_New_Only(flow % p,  grid, 'P')

  ! Allocate memory for mass and volumetric fluxes
  call Face_Mod_Allocate_New_And_Old(flow % m_flux, grid, 'M_FLUX')

  !-----------------------------------------!
  !   Enthalpy conservation (temperature)   !
  !-----------------------------------------!
  if(heat_transfer) then
    call Var_Mod_Allocate_Solution(flow % t, grid, 'T', 'Q')
  end if ! heat_transfer

  allocate(flow % vort (-grid % n_bnd_cells:grid % n_cells)); flow % vort  = 0.
  allocate(flow % shear(-grid % n_bnd_cells:grid % n_cells)); flow % shear = 0.

  allocate(flow % body_fx(-grid % n_bnd_cells:grid % n_cells))
  allocate(flow % body_fy(-grid % n_bnd_cells:grid % n_cells))
  allocate(flow % body_fz(-grid % n_bnd_cells:grid % n_cells))
  flow % body_fx = 0.0
  flow % body_fy = 0.0
  flow % body_fz = 0.0

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
    call Var_Mod_Allocate_Solution(flow % scalar(sc), grid, c_name, q_name)

  end do

  end subroutine
