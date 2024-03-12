#include "../Shared/Assert.h90"
#include "../Shared/Browse.h90"

!==============================================================================!
  module Field_Mod
!----------------------------------[Modules]-----------------------------------!
  use Assert_Mod
  use Bulk_Mod
  use Sparse_Mod
  use Native_Mod
  use Var_mod
  use Gpu_mod
  use Profiler_mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  ! Frequently used aliases (don't know where to put them yet)
  real, contiguous, pointer :: u_n(:), u_o(:)
  real, contiguous, pointer :: v_n(:), v_o(:)
  real, contiguous, pointer :: w_n(:), w_o(:)
  real, contiguous, pointer :: pp_n(:)
  real, contiguous, pointer :: pp_x(:)
  real, contiguous, pointer :: pp_y(:)
  real, contiguous, pointer :: pp_z(:)
  real, contiguous, pointer :: p_n(:)
  real, contiguous, pointer :: p_x(:)
  real, contiguous, pointer :: p_y(:)
  real, contiguous, pointer :: p_z(:)

  integer :: grid_n_cells
  integer :: grid_n_faces
  integer :: grid_n_bnd_cells
  integer :: grid_n_buff_cells

  integer, contiguous, pointer :: grid_faces_c(:,:)
  integer, contiguous, pointer :: grid_cells_c(:,:)
  integer, contiguous, pointer :: grid_cells_f(:,:)
  integer, contiguous, pointer :: grid_cells_n_cells(:)
  integer, contiguous, pointer :: grid_reg_f_face(:)
  integer, contiguous, pointer :: grid_reg_l_face(:)
  real,    contiguous, pointer :: grid_vol(:)
  real,    contiguous, pointer :: grid_xc(:)
  real,    contiguous, pointer :: grid_yc(:)
  real,    contiguous, pointer :: grid_zc(:)
  real,    contiguous, pointer :: grid_dx(:)
  real,    contiguous, pointer :: grid_dy(:)
  real,    contiguous, pointer :: grid_dz(:)
  real,    contiguous, pointer :: grid_sx(:)
  real,    contiguous, pointer :: grid_sy(:)
  real,    contiguous, pointer :: grid_sz(:)
  real,    contiguous, pointer :: grid_s(:)
  real,    contiguous, pointer :: grid_d(:)

  !----------------!
  !                !
  !   Field type   !
  !                !
  !----------------!
  type Field_Type

    type(Grid_Type), pointer :: pnt_grid    !! grid for which it is defined

    !-------------------------!
    !   Physical properties   !
    !-------------------------!

    ! Defined as constants
    real :: capacity     = 1.0     !! [J/kg/K]
    real :: conductivity = 0.01    !! [W/(m K)]
    real :: density      = 1.0     !! [kg / m^3]
    real :: viscosity    = 0.001   !! [kg / (m s)]
    real :: diffusivity  = 1.0e-6  !! [m^2/s]

    !---------------------------------------------------!
    !   Associated with momentum conservation eqution   !
    !---------------------------------------------------!

    ! Velocity components
    type(Var_Type) :: u    !! velocity component [m/s]
    type(Var_Type) :: v    !! velocity component [m/s]
    type(Var_Type) :: w    !! velocity component [m/s]

    ! Pressure-like potential for initial velocity field
    real, allocatable :: potential(:)

    ! Volume flow rate through cell faces
    real, allocatable :: v_flux(:)  !! volume flow rate [m^3/s]

    ! Pressure and pressure correction
    type(Var_Type) :: p    !! pressure            [N/m^2] = [kg/m/s^2]
    type(Var_Type) :: pp   !! pressure correction [N/m^2] = [kg/m/s^2]

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

    ! Time step used in this field
    real :: dt  !! time step in this field [s]

    ! Volume error after pressure correction
    ! (It used to be called mass_err and was a local variable)
    real :: vol_res

    ! Gradient matrices for cells to cells (c2c)
    real, allocatable :: grad_c2c(:,:)  !! gradient matrices [1/m^2]

    ! Some numerical parameters
    real :: blend  !! bleding coefficient for momentum

    contains

      !------------------------!
      !   Core functionality   !
      !------------------------!
      procedure :: Create_Field
      procedure :: Update_Aliases

      !-----------------------------------------!
      !   Procedures for gradient computation   !
      !-----------------------------------------!
      procedure :: Calculate_Grad_Matrix
      procedure :: Grad_Component
      procedure :: Grad_Pressure
      !---------------!
      !   Utilities   !
      !---------------!
      procedure :: Adjust_P_Drops
      procedure :: Alias_Energy
      procedure :: Alias_Momentum
      procedure :: Calculate_Bulk_Fluxes
      procedure :: Volume_Average

  end type

  contains

    !------------------------!
    !   Core functionality   !
    !------------------------!
#   include "Field_Mod/Core/Create_Field.f90"
#   include "Field_Mod/Core/Update_Aliases.f90"

    !-----------------------------------------!
    !   Procedures for gradient computation   !
    !-----------------------------------------!
#   include "Field_Mod/Gradients/Calculate_Grad_Matrix.f90"
#   include "Field_Mod/Gradients/Grad_Component.f90"
#   include "Field_Mod/Gradients/Grad_Pressure.f90"

    !---------------!
    !   Utilities   !
    !---------------!
#   include "Field_Mod/Utilities/Adjust_P_Drops.f90"
#   include "Field_Mod/Utilities/Alias_Energy.f90"
#   include "Field_Mod/Utilities/Alias_Momentum.f90"
#   include "Field_Mod/Utilities/Calculate_Bulk_Fluxes.f90"
#   include "Field_Mod/Utilities/Volume_Average.f90"

  end module
