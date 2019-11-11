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
  use Grid_Mod, only: Grid_Type
  use Bulk_Mod, only: Bulk_Type
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !----------------!
  !   Field type   !
  !----------------!
  type Field_Type

    type(Grid_Type), pointer :: pnt_grid  ! grid for which it is defined

    ! Velocity components
    type(Var_Type) :: u
    type(Var_Type) :: v
    type(Var_Type) :: w

    ! Mass and volumetric flux through cell faces:
    type(Face_Type) :: m_flux
    type(Face_Type) :: v_flux

    ! Pressure and pressure correction
    type(Var_Type) :: p
    type(Var_Type) :: pp

    ! Temperature
    type(Var_Type) :: t

    ! Shear and wall stress are used in a number of turbulence models
    real, allocatable :: shear(:)
    real, allocatable :: vort(:)

    ! Scalars (like chemical species for example)
    integer                     :: n_scalars
    type(Var_Type), allocatable :: scalar(:)

    ! Bulk velocities, pressure drops, etc.
    type(Bulk_Type) :: bulk

    ! Time step used in this field
    real :: dt

    ! Volume expansion coefficient
    real :: beta

  end type

  ! Variables determining if we are dealing with heat transfer and buoyancy
  logical :: heat_transfer
  logical :: buoyancy

  ! Heat flux to the domain (important for periodic case with heat transfer)
  real :: heat_flux, heated_area, heat

  ! Physical properties
  real :: conductivity, diffusivity, capacity
  real, allocatable :: viscosity(:), density(:), dens_face(:)

  ! Angular velocity
  real :: omega_x, omega_y, omega_z, omega

  ! Gravity
  real :: grav_x, grav_y, grav_z

  ! Reference temperature
  real :: t_ref

  contains

  include 'Field_Mod/Allocate.f90'
  include 'Field_Mod/Alias_Energy.f90'
  include 'Field_Mod/Alias_Momentum.f90'
  include 'Field_Mod/U_Tan.f90'

  end module
