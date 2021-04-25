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
  integer       :: sc, nb, nc, nn, nf
  character(VL) :: c_name, q_name
!==============================================================================!

  ! Store the pointer to a grid
  flow % pnt_grid => grid

  ! Take some aliases
  nb = grid % n_bnd_cells
  nc = grid % n_cells
  nn = grid % n_nodes
  nf = grid % n_faces

  !---------------------------------------------!
  !   Allocate memory for physical properties   !
  !---------------------------------------------!
  allocate(flow % density     (-nb:nc));  flow % density(:)      = 0.0
  allocate(flow % viscosity   (-nb:nc));  flow % viscosity(:)    = 0.0
  allocate(flow % capacity    (-nb:nc));  flow % capacity(:)     = 0.0
  allocate(flow % conductivity(-nb:nc));  flow % conductivity(:) = 0.0

  !----------------------------------!
  !   Memory for gradient matrices   !  (are the latter two used at all?)
  !----------------------------------!
  allocate(flow % grad_c2c(6, nc));  flow % grad_c2c(:,:) = 0.0
  allocate(flow % grad_n2c(6, nc));  flow % grad_n2c(:,:) = 0.0
  allocate(flow % grad_c2n(6, nn));  flow % grad_c2n(:,:) = 0.0

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

  ! Allocate memory for volumetric fluxes
  call Face_Mod_Allocate(flow % v_flux, grid, 'V_FL')

  !-----------------------------------------!
  !   Enthalpy conservation (temperature)   !
  !-----------------------------------------!
  if(flow % heat_transfer) then
    call Var_Mod_Allocate_Solution(flow % t, grid, 'T', 'Q')
  end if ! heat_transfer

  allocate(flow % vort (-nb:nc)); flow % vort  = 0.
  allocate(flow % shear(-nb:nc)); flow % shear = 0.

  !----------------------------------------------------------------!
  !   Twelve variables which follow are needed for Rhie and Chow   !
  !----------------------------------------------------------------!
  allocate(flow % fx(nc));           flow % fx = 0.0
  allocate(flow % fy(nc));           flow % fy = 0.0
  allocate(flow % fz(nc));           flow % fz = 0.0

  allocate(flow % cell_fx(-nb:nc));  flow % cell_fx = 0.0
  allocate(flow % cell_fy(-nb:nc));  flow % cell_fy = 0.0
  allocate(flow % cell_fz(-nb:nc));  flow % cell_fz = 0.0

  allocate(flow % face_fx(nf));      flow % face_fx = 0.0
  allocate(flow % face_fy(nf));      flow % face_fy = 0.0
  allocate(flow % face_fz(nf));      flow % face_fz = 0.0

  allocate(flow % u_star(-nb:nc));   flow % u_star      = 0.0
  allocate(flow % v_star(-nb:nc));   flow % v_star      = 0.0
  allocate(flow % w_star(-nb:nc));   flow % w_star      = 0.0
  allocate(flow % v_flux_star(nf));  flow % v_flux_star = 0.0

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
