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

  !----------------!
  !   Field type   !
  !----------------!
  type Field_Type

    type(Grid_Type), pointer :: pnt_grid    !! grid for which it is defined

    !-------------------------!
    !   Physical properties   !
    !-------------------------!

    ! Defined in cell centers
    real, allocatable :: capacity(:)      !! [J/kg/K]
    real, allocatable :: conductivity(:)  !! [W/(m K)]
    real, allocatable :: density(:)       !! [kg/m^3]
    real, allocatable :: viscosity(:)     !! [kg/m/s]
    real              :: diffusivity      !! [m^2/s]

    ! Helping variable for easier forming of conservation equations
    real, allocatable :: ones(:)          !! [1]

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

    ! Used in pressure discretization
    real, allocatable :: v_m(:)    !! cell volume over momentum diagonal

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

    ! Tolerance and maximum iterations for Gauss gradients
    integer :: least_miter
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
      procedure :: Grad_Variable

      !---------------!
      !   Utilities   !
      !---------------!
      procedure :: Adjust_P_Drops
      procedure :: Alias_Energy
      procedure :: Alias_Momentum
      procedure :: Calculate_Bulk_Velocities
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
#   include "Field_Mod/Gradients/Grad_Variable.f90"

    !---------------!
    !   Utilities   !
    !---------------!
#   include "Field_Mod/Utilities/Adjust_P_Drops.f90"
#   include "Field_Mod/Utilities/Alias_Energy.f90"
#   include "Field_Mod/Utilities/Alias_Momentum.f90"
#   include "Field_Mod/Utilities/Calculate_Bulk_Velocities.f90"
#   include "Field_Mod/Utilities/Volume_Average.f90"

  end module
