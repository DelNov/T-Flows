#include "../Shared/Assert.h90"
#include "../Shared/Browse.h90"

!==============================================================================!
  module Field_Mod
!------------------------------------------------------------------------------!
!   Module for basic flow field plus temperature.                              !
!   It is a bit of a mumbo-jumbo at this moment, it will furhter have to       !
!   differentiate into numerical and physica parts.                            !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Face_Mod
  use Bulk_Mod
  use Solver_Mod
  use Time_Mod
  use Numerics_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !----------------!
  !   Field type   !
  !----------------!
  type Field_Type

    type(Matrix_Type), pointer :: pnt_matrix  ! pointer to the matrix
    type(Grid_Type),   pointer :: pnt_grid    ! grid for which it is defined

    !-------------------------!
    !   Physical properties   !
    !-------------------------!

    ! Defined in cell centers
    real, allocatable :: capacity(:)      ! [J/kg/K]
    real, allocatable :: conductivity(:)  ! [W/(m K)]
    real, allocatable :: density(:)       ! [kg/m^3]
    real, allocatable :: viscosity(:)     ! [kg/m/s]
    real              :: diffusivity      ! [m^2/s]

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

    ! Pressure-like potential for initial velocity field
    real, allocatable :: potential(:)

    ! Volumetric flux through cell faces
    type(Face_Type) :: v_flux  ! [m^3/s]

    ! Pressure and pressure correction
    type(Var_Type) :: p    ! [N/m^2] = [kg/m/s^2]
    type(Var_Type) :: pp   ! [N/m^2] = [kg/m/s^2]
    type(Var_Type) :: pot  ! pressure-like potential for initial velocity field

    ! Wall distance used for computation from partial differential equation
    type(Var_Type) :: wall_dist

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

    ! True if it has outlets (needed for a fix in Compute_Pressure)
    !
    ! Update on July 17, 2021: I have some reservations about this part, since
    ! there was another bug fix when computing fluxes in the meanwhile (check:
    ! 90f77a1c8bd4ca05330a4435ed6321782ef00199).  This balancing also caused a
    ! bug when loading backup file (also check "Initialize_Variables",
    ! "Compute_Pressure" and "Backup_Mod/Load and Backup_Mod/Save" procedures)
    !
    ! Update on February 27, 2022: I have also added "has_outflow_boundary"
    ! to be able to tell PETSc if matrix for pressure is singular
    !
    ! Update on June 2, 2022: Unified all outlet boundaries into one
    ! to be able to tell PETSc if matrix for pressure is singular
    logical :: has_pressure

    !-------------------------------------------------!
    !   Associated with energy conservation eqution   !
    !-------------------------------------------------!

    ! Variables determining if we are dealing with heat transfer and buoyancy
    logical :: heat_transfer

    ! Phase change (called mass_transfer to be consistent with heat_transfer)
    ! 0 = no mass transfer, 1 = temperature gradient, 2 = Lee model
    integer :: mass_transfer_model

    ! Temperature
    type(Var_Type) :: t  ! [K]

    ! Heat flux to the domain (important for periodic case with heat transfer)
    real :: heat_flux, heated_area, heat  ! [W/m^2], [m^2], [W]

    ! Reference temperature and volume expansion coefficient (for buoyancy)
    real :: t_ref  ! [K]
    real :: beta   ! [1/K]

    ! Exponential extrapolation of temperature to the walls
    logical :: exp_temp_wall

    ! Scalars (like chemical species for example)
    integer                     :: n_scalars
    type(Var_Type), allocatable :: scalar(:)

    ! Bulk velocities, pressure drops, etc.
    type(Bulk_Type) :: bulk

    !------------------------------------------!
    !   Associated with multiphase situation   !
    !------------------------------------------!
    logical :: with_particles
    logical :: with_interface

    !--------------------------!
    !   Numerical parameters   !
    !--------------------------!

    ! Pressure velocity coupling algorithm
    integer :: p_m_coupling, i_corr, n_piso_corrections
    logical :: inside_piso_loop
    logical :: choi_correction
    logical :: gu_correction
    logical :: rep_vol_balance

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
    integer :: least_miter
    real    :: gauss_iters
    integer :: gauss_calls

    ! Is buoyancy thermally- or density-driven?
    integer :: buoyancy

    ! Gravity must be part of field for conjugate heat transfer models
    real :: grav_x, grav_y, grav_z

    ! Angular velocity
    real :: omega_x, omega_y, omega_z

    ! For volume balance reporting
    integer :: fuvbr

    contains

      !------------------------!
      !   Core functionality   !
      !------------------------!
      procedure :: Create_Field

      !-----------------------------------------!
      !   Procedures for gradient computation   !
      !-----------------------------------------!
      procedure          :: Calculate_Grad_Matrix
      ! procedure          :: Calculate_Grad_Matrix_Cell_By_Cell    ! not used
      ! procedure          :: Calculate_Grad_Matrix_For_Cell        ! not used
      ! procedure          :: Calculate_Grad_Matrix_Faces_To_Cells  ! not used
      ! procedure          :: Calculate_Grad_Matrix_Nodes_To_Cells  ! not used
      ! procedure          :: Calculate_Grad_Matrix_Cells_To_Nodes  ! not used
      procedure          :: Grad
      procedure          :: Grad_Component
      procedure, private :: Grad_Component_No_Refresh
      procedure, private :: Grad_Three_Components_No_Refresh
      ! procedure          :: Grad_Component_Faces_To_Cells         ! not used
      ! procedure          :: Grad_Component_Nodes_To_Cells         ! not used
      ! procedure          :: Grad_Component_Cells_To_Nodes         ! not used
      procedure          :: Grad_Gauss
      procedure, private :: Grad_Gauss_Pressure
      procedure, private :: Grad_Gauss_Variable
      procedure, private :: Grad_Least_Pressure
      procedure, private :: Grad_Least_Variable
      procedure          :: Grad_Pressure
      procedure          :: Grad_Variable

      !----------------------------------!
      !   Procedures for interpolation   !
      !----------------------------------!
      procedure :: Interpolate_Cells_To_Nodes
      ! procedure :: Interpolate_Nodes_To_Cells
      ! procedure :: Interpolate_Nodes_To_Faces
      procedure :: Interpolate_To_Faces_Harmonic
      procedure :: Interpolate_To_Faces_Linear

      !---------------!
      !   Utilities   !
      !---------------!
      procedure :: Alias_Energy
      procedure :: Alias_Momentum
      procedure :: Buoyancy_Forces
      procedure :: Calculate_Courant_In_Cells  ! for post-processing
      procedure :: Calculate_Bulk_Fluxes
      procedure :: Compute_Wall_Distance       ! see: Potential_Initialization
      procedure :: Potential_Initialisation    ! see: Compute_Wall_Distance
      procedure :: Prandtl_Numb
      procedure :: Schmidt_Numb
      procedure :: U_Tan
      procedure :: Report_Vol_Balance
      procedure :: Report_Vol_Balance_Start
      procedure :: Report_Vol_Balance_Stop
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
#   include "Field_Mod/Core/Create_Field.f90"

    !-----------------------------------------!
    !   Procedures for gradient computation   !
    !-----------------------------------------!
#   include "Field_Mod/Gradients/Calculate_Grad_Matrix.f90"
!   include "Field_Mod/Gradients/Calculate_Grad_Matrix_Cell_By_Cell.f90"
!   include "Field_Mod/Gradients/Calculate_Grad_Matrix_For_Cell.f90"
!   include "Field_Mod/Gradients/Calculate_Grad_Matrix_Faces_To_Cells.f90"
!   include "Field_Mod/Gradients/Calculate_Grad_Matrix_Nodes_To_Cells.f90"
!   include "Field_Mod/Gradients/Calculate_Grad_Matrix_Cells_To_Nodes.f90"
#   include "Field_Mod/Gradients/Grad.f90"
#   include "Field_Mod/Gradients/Grad_Component.f90"
#   include "Field_Mod/Gradients/Grad_Component_No_Refresh.f90"
#   include "Field_Mod/Gradients/Grad_Three_Components_No_Refresh.f90"
!   include "Field_Mod/Gradients/Grad_Component_Faces_To_Cells.f90"
!   include "Field_Mod/Gradients/Grad_Component_Nodes_To_Cells.f90"
!   include "Field_Mod/Gradients/Grad_Component_Cells_To_Nodes.f90"
#   include "Field_Mod/Gradients/Grad_Gauss.f90"
#   include "Field_Mod/Gradients/Grad_Gauss_Pressure.f90"
#   include "Field_Mod/Gradients/Grad_Gauss_Variable.f90"
#   include "Field_Mod/Gradients/Grad_Least_Pressure.f90"
#   include "Field_Mod/Gradients/Grad_Least_Variable.f90"
#   include "Field_Mod/Gradients/Grad_Pressure.f90"
#   include "Field_Mod/Gradients/Grad_Variable.f90"

    !----------------------------------!
    !   Procedures for interpolation   !
    !----------------------------------!
#   include "Field_Mod/Interpolations/Interpolate_Cells_To_Nodes.f90"
!   include "Field_Mod/Interpolations/Interpolate_Nodes_To_Cells.f90"
!   include "Field_Mod/Interpolations/Interpolate_Nodes_To_Faces.f90"
#   include "Field_Mod/Interpolations/Interpolate_To_Faces_Harmonic.f90"
#   include "Field_Mod/Interpolations/Interpolate_To_Faces_Linear.f90"

    !---------------!
    !   Utilities   !
    !---------------!
#   include "Field_Mod/Utilities/Alias_Energy.f90"
#   include "Field_Mod/Utilities/Alias_Momentum.f90"
#   include "Field_Mod/Utilities/Buoyancy_Forces.f90"
#   include "Field_Mod/Utilities/Calculate_Courant_In_Cells.f90"
#   include "Field_Mod/Utilities/Calculate_Bulk_Fluxes.f90"
#   include "Field_Mod/Utilities/Calculate_Shear_And_Vorticity.f90"
#   include "Field_Mod/Utilities/Potential_Initialization.f90"
#   include "Field_Mod/Utilities/Prandtl_Numb.f90"
#   include "Field_Mod/Utilities/Schmidt_Numb.f90"
#   include "Field_Mod/Utilities/U_Tan.f90"
#   include "Field_Mod/Utilities/Compute_Wall_Distance.f90"
#   include "Field_Mod/Utilities/Report_Vol_Balance.f90"
#   include "Field_Mod/Utilities/Report_Vol_Balance_Start.f90"
#   include "Field_Mod/Utilities/Report_Vol_Balance_Stop.f90"
#   include "Field_Mod/Utilities/Volume_Average.f90"

  end module
