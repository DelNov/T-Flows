!==============================================================================!
  subroutine Create_Field(Flow, A)
!------------------------------------------------------------------------------!
!   Allocates memory for the entire field.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Field_Type)         :: Flow
  type(Matrix_Type), target :: A
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
  integer                  :: sc, nb, nc, nn, nf
  character(VL)            :: c_name, q_name
!==============================================================================!

  ! Store the pointer to a Grid
  Flow % pnt_matrix => A
  Flow % pnt_grid   => A % pnt_grid
  Grid              => A % pnt_grid

  ! Take some aliases
  nb = Grid % n_bnd_cells
  nc = Grid % n_cells
  nn = Grid % n_nodes
  nf = Grid % n_faces

  !---------------------------------------------!
  !   Allocate memory for physical properties   !
  !---------------------------------------------!
  allocate(Flow % density     (-nb:nc));  Flow % density(:)      = 0.0
  allocate(Flow % viscosity   (-nb:nc));  Flow % viscosity(:)    = 0.0
  allocate(Flow % capacity    (-nb:nc));  Flow % capacity(:)     = 0.0
  allocate(Flow % conductivity(-nb:nc));  Flow % conductivity(:) = 0.0

  !----------------------------------!
  !   Memory for gradient matrices   !  (are the latter two used at all?)
  !----------------------------------!
  allocate(Flow % grad_c2c(6, nc));  Flow % grad_c2c(:,:) = 0.0
  allocate(Flow % grad_f2c(6, nc));  Flow % grad_f2c(:,:) = 0.0
  allocate(Flow % grad_n2c(6, nc));  Flow % grad_n2c(:,:) = 0.0
  allocate(Flow % grad_c2n(6, nn));  Flow % grad_c2n(:,:) = 0.0

  !----------------------------!
  !   Navier-Stokes equation   !
  !----------------------------!

  ! Allocate memory for velocity components
  call Var_Mod_Create_Solution(Flow % u, A, 'U', '')
  call Var_Mod_Create_Solution(Flow % v, A, 'V', '')
  call Var_Mod_Create_Solution(Flow % w, A, 'W', '')

  ! Potential for initial velocity computation
  call Var_Mod_Create_Solution(Flow % pot, A, 'POT', '')

  ! For computation of wall distance
  call Var_Mod_Create_Solution(Flow % wall_dist, A, 'WD', '')

  ! Allocate memory for pressure correction and pressure
  call Var_Mod_Create_Solution(Flow % pp, A, 'PP', '')
  call Var_Mod_Create_New_Only(Flow % p,  Grid, 'P')

  ! Allocate memory for volumetric fluxes
  call Face_Mod_Allocate(Flow % v_flux, Grid, 'V_FL')

  !-----------------------------------------!
  !   Enthalpy conservation (temperature)   !
  !-----------------------------------------!
  if(Flow % heat_transfer) then
    call Var_Mod_Create_Solution(Flow % t, A, 'T', 'Q')
  end if ! heat_transfer

  allocate(Flow % vort (-nb:nc)); Flow % vort  = 0.
  allocate(Flow % shear(-nb:nc)); Flow % shear = 0.

  !--------------------------------------------------------------!
  !   Nine variables which follow are needed for Rhie and Chow   !
  !--------------------------------------------------------------!
  allocate(Flow % fx(nc));           Flow % fx = 0.0
  allocate(Flow % fy(nc));           Flow % fy = 0.0
  allocate(Flow % fz(nc));           Flow % fz = 0.0

  allocate(Flow % cell_fx(-nb:nc));  Flow % cell_fx = 0.0
  allocate(Flow % cell_fy(-nb:nc));  Flow % cell_fy = 0.0
  allocate(Flow % cell_fz(-nb:nc));  Flow % cell_fz = 0.0

  allocate(Flow % face_fx(nf));      Flow % face_fx = 0.0
  allocate(Flow % face_fy(nf));      Flow % face_fy = 0.0
  allocate(Flow % face_fz(nf));      Flow % face_fz = 0.0

  !--------------------------------------!
  !   Allocate memory for user scalars   !
  !--------------------------------------!
  allocate(Flow % scalar(Flow % n_scalars))

  !-------------------------------------!
  !   Browse through all user scalars   !
  !-------------------------------------!
  do sc = 1, Flow % n_scalars

    ! Set variable name
    c_name = 'C_00'
    q_name = 'Q_00'
    write(c_name(3:4),'(i2.2)') sc
    write(q_name(3:4),'(i2.2)') sc

    ! Allocate memory for passive scalar
    call Var_Mod_Create_Solution(Flow % scalar(sc), A, c_name, q_name)

  end do

  !--------------------------------------------------------!
  !   Initialize counters for Gauss gradient computation   !
  !--------------------------------------------------------!
  Flow % gauss_iters = 0.0
  Flow % gauss_calls = 0

  end subroutine
