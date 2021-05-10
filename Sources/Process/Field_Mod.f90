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
  use Numerics_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !----------------!
  !   Field type   !
  !----------------!
  type Field_Type

    type(Grid_Type), pointer :: pnt_grid  ! grid for which it is defined

    !-------------------------!
    !   Physical properties   !
    !-------------------------!

    ! Defined in cell centers
    real, allocatable :: capacity(:)      ! [J/kg/K]
    real, allocatable :: conductivity(:)  ! [W/(m K)]
    real, allocatable :: density(:)       ! [kg/m^3]
    real, allocatable :: viscosity(:)     ! [kg/m/s]

    ! Defined globally
    real :: diffusivity      ! [m^2/s]
    real :: latent_heat      ! [J/kg]
    real :: sat_temperature  ! [K]

    !---------------------------------------------------!
    !   Associated with momentum conservation eqution   !
    !---------------------------------------------------!

    ! Velocity components
    type(Var_Type) :: u    ! [m/s]
    type(Var_Type) :: v    ! [m/s]
    type(Var_Type) :: w    ! [m/s]

    ! Shear and wall stress are used in a number of turbulence models
    real, allocatable :: shear(:)  ! [1/s]
    real, allocatable :: vort(:)   ! [1/s]

    ! Volumetric flux through cell faces
    type(Face_Type) :: v_flux  ! [m^3/s]

    ! Pressure and pressure correction
    type(Var_Type) :: p    ! [N/m^2] = [kg/m/s^2]
    type(Var_Type) :: pp   ! [N/m^2] = [kg/m/s^2]
    type(Var_Type) :: pot  ! pressure-like potential for initial velocity field

    ! Internal forces on the fluid.
    ! These includes forces due to discretization (cross diffusion terms),
    ! boundary conditions (diffusive or advective), under-relaxation
    real, allocatable :: fx(:)
    real, allocatable :: fy(:)
    real, allocatable :: fz(:)

    ! External "body" forces on the fluid.
    ! These includes body forces such as: buoyancy, gravity, surface tension
    ! These need special treatment according to Mencinger and Zun (2007) for
    ! proper linking with Rhie and Chow algorithm.  Therefore, two variants
    ! are included: face based forces and cell based forces.  The cell based
    ! forces must be derived from face based forces and properly linked.
    real, allocatable :: cell_fx(:)  ! [N/m^3]
    real, allocatable :: cell_fy(:)  ! [N/m^3]
    real, allocatable :: cell_fz(:)  ! [N/m^3]
    real, allocatable :: face_fx(:)  ! [N/m^3]
    real, allocatable :: face_fy(:)  ! [N/m^3]
    real, allocatable :: face_fz(:)  ! [N/m^3]

    ! Reference density (for buoyancy)
    real :: dens_ref

    !-------------------------------------------------!
    !   Associated with energy conservation eqution   !
    !-------------------------------------------------!

    ! Variables determining if we are dealing with heat transfer and buoyancy
    logical :: heat_transfer

    ! Phase change (called mass_transfer to be consistent with heat_transfer)
    logical :: mass_transfer

    ! Temperature
    type(Var_Type) :: t  ! [K]

    ! Heat flux to the domain (important for periodic case with heat transfer)
    real :: heat_flux, heated_area, heat  ! [W/m^2], [m^2], [W]

    ! Reference temperature and volume expansion coefficient (for buoyancy)
    real :: t_ref  ! [K]
    real :: beta   ! [1/K]

    ! Scalars (like chemical species for example)
    integer                     :: n_scalars
    type(Var_Type), allocatable :: scalar(:)

    ! Bulk velocities, pressure drops, etc.
    type(Bulk_Type) :: bulk

    !--------------------------!
    !   Numerical parameters   !
    !--------------------------!

    ! Pressure velocity coupling algorithm
    integer :: p_m_coupling, i_corr, n_piso_corrections
    logical :: piso_status
    logical :: choi_correction
    logical :: gu_correction

    ! Maximum CFL and Pe numbers
    real :: cfl_max, pe_max

    ! Time step used in this field
    real :: dt  ! [s]

    ! Volume error after pressure correction
    ! (It used to be called mass_err and was a local variable)
    real :: vol_res

    ! Gradient matrices for:
    ! - cells to cells (c2c)
    ! - nodes to cells (n2c), and
    ! - cells to nodes (c2n)
    real, allocatable :: grad_c2c(:,:)
    real, allocatable :: grad_f2c(:,:)
    real, allocatable :: grad_n2c(:,:)
    real, allocatable :: grad_c2n(:,:)

    ! Tolerance and maximum iterations for Gauss gradients
    real    :: gauss_tol
    integer :: gauss_miter

    ! Is buoyancy thermally- or density-driven?
    integer :: buoyancy

  end type

  ! Angular velocity
  real :: omega_x, omega_y, omega_z, omega

  ! Gravity
  real, target :: grav_x, grav_y, grav_z

  ! Parameters for type of buoyancy
  integer, parameter :: NO_BUOYANCY      = 60013
  integer, parameter :: DENSITY_DRIVEN   = 60017
  integer, parameter :: THERMALLY_DRIVEN = 60029

  contains

  include 'Field_Mod/Allocate.f90'
  include 'Field_Mod/Alias_Energy.f90'
  include 'Field_Mod/Alias_Momentum.f90'
  include 'Field_Mod/Buoyancy_Forces.f90'
  include 'Field_Mod/Calculate_Fluxes.f90'
  include 'Field_Mod/Calculate_Grad_Matrix.f90'
  include 'Field_Mod/Calculate_Grad_Matrix_Cell_By_Cell.f90'
  include 'Field_Mod/Calculate_Grad_Matrix_For_Cell.f90'
  include 'Field_Mod/Calculate_Grad_Matrix_Faces_To_Cells.f90'
  include 'Field_Mod/Calculate_Grad_Matrix_Nodes_To_Cells.f90'
  include 'Field_Mod/Calculate_Grad_Matrix_Cells_To_Nodes.f90'
  include 'Field_Mod/Grad.f90'
  include 'Field_Mod/Grad_Component.f90'
  include 'Field_Mod/Grad_Component_No_Refresh.f90'
  include 'Field_Mod/Grad_Component_Faces_To_Cells.f90'
  include 'Field_Mod/Grad_Component_Nodes_To_Cells.f90'
  include 'Field_Mod/Grad_Component_Cells_To_Nodes.f90'
  include 'Field_Mod/Grad_Gauss.f90'
  include 'Field_Mod/Grad_Gauss_Pressure.f90'
  include 'Field_Mod/Grad_Gauss_Variable.f90'
  include 'Field_Mod/Grad_Least_Pressure.f90'
  include 'Field_Mod/Grad_Least_Pressure_Correction.f90'
  include 'Field_Mod/Grad_Least_Variable.f90'
  include 'Field_Mod/Grad_Pressure.f90'
  include 'Field_Mod/Grad_Pressure_Correction.f90'
  include 'Field_Mod/Grad_Variable.f90'
  include 'Field_Mod/Interpolate_Cells_To_Nodes.f90'
  include 'Field_Mod/Interpolate_Nodes_To_Cells.f90'
  include 'Field_Mod/Interpolate_Nodes_To_Faces.f90'
  include 'Field_Mod/Interpolate_To_Faces_Harmonic.f90'
  include 'Field_Mod/Interpolate_To_Faces_Linear.f90'
  include 'Field_Mod/Potential_Initialization.f90'
  include 'Field_Mod/Prandtl_Number.f90'
  include 'Field_Mod/Schmidt_Number.f90'
  include 'Field_Mod/U_Tan.f90'

  end module
