!==============================================================================!
  module Field_Mod
!------------------------------------------------------------------------------!
!   Module for basic flow field plus temperature.                              !
!   It is a bit of a mumbo-jumbo at this moment, it will furhter have to       !
!   differentiate into numerical and physica parts.                            !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Var_Mod
  use Face_Mod
  use Grid_Mod
  use Bulk_Mod
  use Solver_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !----------------!
  !   Field type   !
  !----------------!
  type Field_Type

    type(Grid_Type), pointer :: pnt_grid  ! grid for which it is defined

    ! Pressure velocity coupling algorithm
    integer :: p_m_coupling, i_corr, n_piso_corrections
    logical :: piso_status

    ! Flag for temporal correction for pressure
    logical :: temp_corr

    ! Time step and sub-relaxation coefficient for pressure correction equation
    ! (If temp_corr is .true., the following two variables are needed as well)
    real :: dt_corr, u_rel_corr

    ! Physical properties
    real, allocatable :: capacity(:)      ! [J/kg/K]
    real, allocatable :: conductivity(:)  ! [W/(m K)]
    real, allocatable :: density(:)       ! [kg/m^3]
    real, allocatable :: viscosity(:)     ! [kg/m/s]
    real              :: diffusivity      ! [m^2/s]
    real              :: latent_heat      ! [J/kg]
    real              :: sat_temperature  ! [K]

    ! Velocity components
    type(Var_Type) :: u    ! [m/s]
    type(Var_Type) :: v    ! [m/s]
    type(Var_Type) :: w    ! [m/s]
    type(Var_Type) :: pot  ! potential for initial velocity field

    ! Volumetric flux through cell faces
    type(Face_Type) :: v_flux  ! [m^3/s]

    ! Pressure and pressure correction
    type(Var_Type) :: p   ! [N/m^2] = [kg/m/s^2]
    type(Var_Type) :: pp  ! [N/m^2] = [kg/m/s^2]

    ! Temperature
    type(Var_Type) :: t  ! [K]

    ! Shear and wall stress are used in a number of turbulence models
    real, allocatable :: shear(:)
    real, allocatable :: vort(:)

    ! Scalars (like chemical species for example)
    integer                     :: n_scalars
    type(Var_Type), allocatable :: scalar(:)

    ! Bulk velocities, pressure drops, etc.
    type(Bulk_Type) :: bulk

    ! Maximum CFL and Pe numbers
    real :: cfl_max, pe_max

    ! Time step used in this field
    real :: dt  ! [s]

    ! Volume expansion coefficient
    real :: beta

    ! Heat flux to the domain (important for periodic case with heat transfer)
    real :: heat_flux, heated_area, heat  ! [W/m^2], [m^2], [W]

    !---------------------------------!
    !   Gradient matrices for:        !
    !   - cells to cells (c2c),       !
    !   - nodes to cells (n2c), and   !
    !   - cells to nodes (c2n)        !
    !---------------------------------!
    real, allocatable :: grad_c2c(:,:)
    real, allocatable :: grad_n2c(:,:)
    real, allocatable :: grad_c2n(:,:)

    ! Body force
    real, allocatable :: body_fx(:)
    real, allocatable :: body_fy(:)
    real, allocatable :: body_fz(:)

    ! Reference temperature and density
    real :: t_ref
    real :: dens_ref

    ! Volume error after pressure correction
    ! (It used to be called mass_err and was a local variable)
    real :: vol_res

  end type

  ! Variables determining if we are dealing with heat transfer and buoyancy
  logical :: heat_transfer

  ! Phase change (called mass_transfer to be consistent with heat_transfer)
  logical :: mass_transfer

  logical :: buoyancy

  ! Angular velocity
  real :: omega_x, omega_y, omega_z, omega

  ! Gravity
  real :: grav_x, grav_y, grav_z

  contains

  include 'Field_Mod/Allocate.f90'
  include 'Field_Mod/Allocate_Grad_Matrix.f90'
  include 'Field_Mod/Alias_Energy.f90'
  include 'Field_Mod/Alias_Momentum.f90'
  include 'Field_Mod/Body_Forces.f90'
  include 'Field_Mod/Calculate_Fluxes.f90'
  include 'Field_Mod/Calculate_Grad_Matrix.f90'
  include 'Field_Mod/Calculate_Grad_Matrix_Cell_By_Cell.f90'
  include 'Field_Mod/Calculate_Grad_Matrix_For_Cell.f90'
  include 'Field_Mod/Calculate_Grad_Matrix_Nodes_To_Cells.f90'
  include 'Field_Mod/Calculate_Grad_Matrix_Cells_To_Nodes.f90'
  include 'Field_Mod/Correct_Fluxes_With_Body_Forces.f90'
  include 'Field_Mod/Grad.f90'
  include 'Field_Mod/Grad_Component.f90'
  include 'Field_Mod/Grad_Component_No_Refresh.f90'
  include 'Field_Mod/Grad_Component_Nodes_To_Cells.f90'
  include 'Field_Mod/Grad_Component_Cells_To_Nodes.f90'
  include 'Field_Mod/Grad_Pressure.f90'
  include 'Field_Mod/Grad_Pressure_Correction.f90'
  include 'Field_Mod/Grad_Variable.f90'
  include 'Field_Mod/Interpolate_Cells_To_Nodes.f90'
  include 'Field_Mod/Interpolate_Nodes_To_Cells.f90'
  include 'Field_Mod/Interpolate_Var_To_Face.f90'
  include 'Field_Mod/Interpolate_To_Face.f90'
  include 'Field_Mod/Potential_Initialization.f90'
  include 'Field_Mod/Prandtl_Number.f90'
  include 'Field_Mod/Schmidt_Number.f90'
  include 'Field_Mod/U_Tan.f90'

  end module
