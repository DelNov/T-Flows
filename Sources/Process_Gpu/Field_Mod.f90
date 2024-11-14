#include "../Shared/Assert.h90"
#include "../Shared/Browse.h90"
#include "../Shared/Macros.h90"
#include "../Shared/Unused.h90"

!==============================================================================!
  module Field_Mod
!----------------------------------[Modules]-----------------------------------!
  use Assert_Mod
  use Face_Mod
  use Bulk_Mod
  use Gpu_Pointers_Mod
  use Sparse_Mod
  use Native_Mod
  use Var_mod
  use Numerics_Mod
  use Profiler_mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !----------------!
  !   Field type   !
  !----------------!
  type Field_Type

    !-------------------------!
    !   Physical properties   !
    !-------------------------!

    ! Defined in cell centers
    real, allocatable :: capacity(:)      !! [J/kg/K]
    real, allocatable :: conductivity(:)  !! [W/(m K)]
    real, allocatable :: density(:)       !! [kg/m^3]
    real, allocatable :: viscosity(:)     !! [kg/m/s]
    real              :: diffusivity      !! [m^2/s]

    !---------------------------------------------------!
    !   Associated with momentum conservation eqution   !
    !---------------------------------------------------!

    ! Velocity components
    type(Var_Type) :: u    !! velocity component [m/s]
    type(Var_Type) :: v    !! velocity component [m/s]
    type(Var_Type) :: w    !! velocity component [m/s]

    ! Shear and wall stress are used in a number of turbulence models
    real, allocatable :: shear(:)  !! shear [1/s]
    real, allocatable :: vort(:)   !! vorticity [1/s]

    ! Pressure-like potential for initial velocity field
    real, allocatable :: potential(:)

    ! Volume flow rate through cell faces
    type(Face_Type) :: v_flux  !! volume flow rate [m^3/s]

    ! Pressure and pressure correction
    type(Var_Type) :: p    !! pressure            [N/m^2] = [kg/m/s^2]
    type(Var_Type) :: pp   !! pressure correction [N/m^2] = [kg/m/s^2]

    ! Used in pressure discretization
    real, allocatable :: v_m(:)    !! cell volume over momentum diagonal

    ! Reference density (for buoyancy)
    real :: dens_ref  !! reference density, used to compute buoyancy terms


    !-------------------------------------------------!
    !   Associated with energy conservation eqution   !
    !-------------------------------------------------!

    ! Variables determining if we are dealing with heat transfer and buoyancy
    logical :: heat_transfer  !! true if heat transfer is invloved

    ! Temperature
    type(Var_Type) :: t  !! temperature [K]

    type(Native_Type) :: Nat

    ! Reference temperature and volume expansion coefficient (for buoyancy)
    real :: t_ref  !! reference temperature used for buoyancy terms [K]
    real :: beta   !! volume expansion coefficient [1/K]

    ! Scalars (like chemical species for example)
    integer                     :: n_scalars  !! number of (passive) scalars
    type(Var_Type), allocatable :: scalar(:)  !! storage for scalars

    ! Bulk velocities, pressure drops, etc.
    type(Bulk_Type) :: bulk  !! holder of volume flow rates through domain and,
                             !! associated with that, bulk velocities

    !--------------------------!
    !   Numerical parameters   !
    !--------------------------!

    ! Pressure velocity coupling algorithm
    logical :: rep_vol_balance

    ! Maximum CFL and Pe numbers
    real :: cfl_max, pe_max

    ! Time step used in this field
    real :: dt  !! time step in this field [s]

    ! Volume error after pressure correction
    ! (It used to be called mass_err and was a local variable)
    real :: vol_res

    ! Gradient matrices for cells to cells (c2c)
    real, allocatable :: grad_c2c(:,:)  !! gradient matrices [1/m^2]
    real, allocatable :: phi_x(:)       !! gradient of phi in x direction
    real, allocatable :: phi_y(:)       !! gradient of phi in y direction
    real, allocatable :: phi_z(:)       !! gradient of phi in z direction
    character(VL) :: stores_gradients_of  !! name of the variable whose
                                          !! gradients are now computed

    ! Tolerance and maximum iterations for Gauss gradients
    integer :: least_miter

    ! Is buoyancy thermally- or density-driven?
    integer :: buoyancy  !! indicates if buoyancy is modeled by Boussinesq
                         !! approach (thermally-driven) or by variations in
                         !! density (density-driven)

    ! Gravity must be part of field for conjugate heat transfer models
    real :: grav_x, grav_y, grav_z

    contains

      !------------------------!
      !   Core functionality   !
      !------------------------!
      procedure :: Add_Advection_Term
      procedure :: Add_Cross_Diffusion_Term
      procedure :: Add_Inertial_Term
      procedure :: Create_Field

      !-----------------------------------------!
      !   Procedures for gradient computation   !
      !-----------------------------------------!
      procedure :: Calculate_Grad_Matrix
      procedure :: Grad_Component
      procedure :: Grad_Pressure
      procedure :: Grad_Variable

      !---------------!
      !   Utilities   !
      !---------------!
      procedure :: Adjust_P_Drops
      procedure :: Alias_Energy
      procedure :: Alias_Momentum
      procedure :: Buoyancy_Forces
      procedure :: Calculate_Bulk_Velocities
      procedure :: Calculate_Shear_And_Vorticity
      procedure :: Volume_Average

  end type

  ! Parameters for type of buoyancy
  integer, parameter :: NO_BUOYANCY      = 60013
  integer, parameter :: DENSITY_DRIVEN   = 60017
  integer, parameter :: THERMALLY_DRIVEN = 60029

  contains

    !------------------------!
    !   Core functionality   !
    !------------------------!
#   include "Field_Mod/Core/Add_Advection_Term.f90"
#   include "Field_Mod/Core/Add_Cross_Diffusion_Term.f90"
#   include "Field_Mod/Core/Add_Inertial_Term.f90"
#   include "Field_Mod/Core/Create_Field.f90"

    !-----------------------------------------!
    !   Procedures for gradient computation   !
    !-----------------------------------------!
#   include "Field_Mod/Gradients/Calculate_Grad_Matrix.f90"
#   include "Field_Mod/Gradients/Grad_Component.f90"
#   include "Field_Mod/Gradients/Grad_Pressure.f90"
#   include "Field_Mod/Gradients/Grad_Variable.f90"

    !---------------!
    !   Utilities   !
    !---------------!
#   include "Field_Mod/Utilities/Adjust_P_Drops.f90"
#   include "Field_Mod/Utilities/Alias_Energy.f90"
#   include "Field_Mod/Utilities/Alias_Momentum.f90"
#   include "Field_Mod/Utilities/Buoyancy_Forces.f90"
#   include "Field_Mod/Utilities/Calculate_Bulk_Velocities.f90"
#   include "Field_Mod/Utilities/Calculate_Shear_And_Vorticity.f90"
#   include "Field_Mod/Utilities/Volume_Average.f90"

  end module
