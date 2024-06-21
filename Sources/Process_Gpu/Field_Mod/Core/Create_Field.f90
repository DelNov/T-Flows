!==============================================================================!
  subroutine Create_Field(Flow, Grid)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Field_Type), target :: Flow  !! parent flow object
  type(Grid_Type),   target :: Grid  !! grid on which the flow is defined
!-----------------------------------[Locals]-----------------------------------!
  integer :: nb, nc, ns
!==============================================================================!

  ! Take some aliases
  nb = Grid % n_bnd_cells
  nc = Grid % n_cells
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

  !--------------------------------------------------------------!
  !   Create native solvers (matrices A, M, right hand side b,   !
  !   helping vectors for CG method such as p, q, r and d_inv)   !
  !--------------------------------------------------------------!
  O_Print '(a)', ' # Creating field''s native solver'
  call Flow % Nat % Create_Native(Grid)

  !----------------------------------!
  !   Memory for gradient matrices   !
  !----------------------------------!
  allocate(Flow % grad_c2c(6, nc));  Flow % grad_c2c(:,:) = 0.0
  allocate(Flow % phi_x(-nb:nc));    Flow % phi_x(:) = 0.0
  allocate(Flow % phi_y(-nb:nc));    Flow % phi_y(:) = 0.0
  allocate(Flow % phi_z(-nb:nc));    Flow % phi_z(:) = 0.0

  !----------------------------!
  !   Navier-Stokes equation   !
  !----------------------------!

  ! Allocate memory for velocity components
  call Var_Mod_Create_Variable(Flow % u, Grid, 'U', '')
  call Var_Mod_Create_Variable(Flow % v, Grid, 'V', '')
  call Var_Mod_Create_Variable(Flow % w, Grid, 'W', '')

  ! Allocate memory for pressure correction and pressure
  call Var_Mod_Create_Variable(Flow % pp, Grid, 'PP', '')
  call Var_Mod_Create_Variable(Flow % p,  Grid, 'P',  '')

  ! Helping array to discretize pressure Poisson equation
  allocate(Flow % v_m(nc));  Flow % v_m(:) = 0.0

  ! Allocate memory for volumetric fluxes
  call Face_Mod_Allocate(Flow % v_flux, Grid, 'V_FL')

  !-----------------------------------------!
  !   Enthalpy conservation (temperature)   !
  !-----------------------------------------!
  if(Flow % heat_transfer) then
    call Var_Mod_Create_Variable(Flow % t, Grid, 'T', 'Q')
  end if ! heat_transfer

  !--------------------------------------------------------------!
  !   Create native solvers (matrices A, M, right hand side b,   !
  !   helping vectors for CG method such as p, q, r and d_inv)   !
  !--------------------------------------------------------------!
  O_Print '(a)', ' # Creating field''s native solver'
  call Flow % Nat % Create_Native(Grid)

  end subroutine
