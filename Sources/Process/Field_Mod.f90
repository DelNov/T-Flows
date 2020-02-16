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
  use Comm_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !----------------!
  !   Field type   !
  !----------------!
  type Field_Type

    type(Grid_Type), pointer :: pnt_grid  ! grid for which it is defined

    ! Physical properties
    real, allocatable :: capacity(:)      ! [J/kg/K]
    real, allocatable :: conductivity(:)  ! [W/(m K)]
    real :: diffusivity                   ! [m^2/s]
    real, allocatable :: density(:)       ! [kg/m^3]
    real, allocatable :: viscosity(:)     ! [kg/m/s]
    real, allocatable :: density_f(:)     ! [kg/m^3]

    ! Velocity components
    type(Var_Type) :: u  ! [m/s]
    type(Var_Type) :: v  ! [m/s]
    type(Var_Type) :: w  ! [m/s]

    ! Mass and volumetric flux through cell faces:
    type(Face_Type) :: m_flux  ! [kg/s]

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

    ! Time step used in this field
    real :: dt  ! [s]

    ! Volume expansion coefficient
    real :: beta

    ! Heat flux to the domain (important for periodic case with heat transfer)
    real :: heat_flux, heated_area, heat  ! [W/m^2], [m^2], [W]

    !-------------------------------------!
    !   Gradient matrices for all cells   !
    !-------------------------------------!
    real, allocatable :: grad(:,:)

  end type

  ! Variables determining if we are dealing with heat transfer and buoyancy
  logical :: heat_transfer
  logical :: buoyancy

  ! Angular velocity
  real :: omega_x, omega_y, omega_z, omega

  ! Gravity
  real :: grav_x, grav_y, grav_z

  ! Reference temperature
  real :: t_ref

  contains

  include 'Field_Mod/Allocate.f90'
  include 'Field_Mod/Allocate_Grad_Matrix.f90'
  include 'Field_Mod/Alias_Energy.f90'
  include 'Field_Mod/Alias_Momentum.f90'
  include 'Field_Mod/Calculate_Fluxes.f90'
  include 'Field_Mod/Calculate_Grad_Matrix.f90'
  include 'Field_Mod/Grad_Component.f90'
  include 'Field_Mod/Grad_Pressure.f90'
  include 'Field_Mod/Grad_Pressure_Correction.f90'
  include 'Field_Mod/Grad_Variable.f90'
  include 'Field_Mod/Prandtl_Number.f90'
  include 'Field_Mod/U_Tan.f90'

  end module
