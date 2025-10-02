!==============================================================================!
  subroutine Create_Field(Flow, A)
!------------------------------------------------------------------------------!
!>  Allocates memory for the entire field in Process.
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * Establishes a pointer to the grid associated with the flow field.        !
!   * Allocates memory for physical properties like density, viscosity,        !
!     capacity, and conductivity.                                              !
!   * Allocates memory for gradient matrices needed for calculations of        !
!     gradients in the flow field.                                             !
!   * Sets up the variables for the solution of Navier-Stokes equations,       !
!     including velocity components, potential for initial velocity            !
!     computation, pressure, and pressure correction.                          !
!   * Handles specific allocations for enthalpy conservation (temperature)     !
!     if heat transfer is involved in the simulation.                          !
!   * Initializes variables related to the Rhie and Chow interpolation method, !
!     such as forces on fluid cells and faces.                                 !
!   * Allocates memory for passive scalars if they are used in the simulation. !
!   * Initializes counters for calculations involving Gauss gradient method.   !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Field_Type), target :: Flow  !! parent flow object
  type(Matrix_Type), target :: A     !! matrix object used with this field
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
  type(Var_Type),  pointer :: u, v, w
  integer                  :: sc, nb, nc, nn, ns
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
  ns = Grid % n_faces

  !---------------------------------------------!
  !   Allocate memory for physical properties   !
  !---------------------------------------------!
  allocate(Flow % density  (-nb:nc));  Flow % density(:)   = 0.0
  allocate(Flow % viscosity(-nb:nc));  Flow % viscosity(:) = 0.0
  if(Flow % heat_transfer) then
    allocate(Flow % capacity    (-nb:nc));  Flow % capacity(:)     = 0.0
    allocate(Flow % conductivity(-nb:nc));  Flow % conductivity(:) = 0.0
  end if
  if(Flow % n_scalars .gt. 0) then
    allocate(Flow % diffusivity (-nb:nc));  Flow % diffusivity(:)  = 0.0
  end if

  !----------------------------------!
  !   Memory for gradient matrices   !
  !----------------------------------!
  allocate(Flow % grad_c2c(6, nc));  Flow % grad_c2c(:,:) = 0.0

  !----------------------------!
  !   Navier-Stokes equation   !
  !----------------------------!
  u => Flow % u
  v => Flow % v
  w => Flow % w

  ! Allocate memory for velocity components
  call Var_Mod_Create_Solution(u, A, 'U', '')
  call Var_Mod_Create_Solution(v, A, 'V', '', reuse_pet = u % pet_rank)
  call Var_Mod_Create_Solution(w, A, 'W', '', reuse_pet = u % pet_rank)

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

  allocate(Flow % potential(-nb:nc)); Flow % potential = 0.

  !--------------------------------------------------------------!
  !   Nine variables which follow are needed for Rhie and Chow   !
  !--------------------------------------------------------------!
  allocate(Flow % fx(nc));           Flow % fx = 0.0
  allocate(Flow % fy(nc));           Flow % fy = 0.0
  allocate(Flow % fz(nc));           Flow % fz = 0.0

  allocate(Flow % cell_fx(-nb:nc));  Flow % cell_fx = 0.0
  allocate(Flow % cell_fy(-nb:nc));  Flow % cell_fy = 0.0
  allocate(Flow % cell_fz(-nb:nc));  Flow % cell_fz = 0.0

  allocate(Flow % face_fx(ns));      Flow % face_fx = 0.0
  allocate(Flow % face_fy(ns));      Flow % face_fy = 0.0
  allocate(Flow % face_fz(ns));      Flow % face_fz = 0.0

  !-----------------------!
  !   Scalars Transport   !
  !-----------------------!
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
