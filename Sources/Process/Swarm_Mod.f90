!==============================================================================!
  module Swarm_Mod
!------------------------------------------------------------------------------!
!   Module for Lagrangian particle tracking                                    !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Particle_Mod
  use Turb_Mod
  use Vof_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  ! Parameters describing turbulence model choice
  ! (Prime numbers starting from 20000)
  integer, parameter :: BROWNIAN_FUKAGATA    = 20011
  integer, parameter :: DISCRETE_RANDOM_WALK = 20021

  !----------------!
  !   Swarm type   !
  !----------------!
  type Swarm_Type

    type(Grid_Type),  pointer :: pnt_grid  ! grid for which it is defined
    type(Field_Type), pointer :: pnt_flow  ! flow field for which it is defined
    type(Turb_Type),  pointer :: pnt_turb  ! turb field for which it is defined
    type(Vof_Type),   pointer :: pnt_vof   ! vof, for three-phase situations

    integer                          :: n_particles = 0
    type(Particle_Type), allocatable :: Particle(:)

    ! Density of this swarm
    real :: density

    ! (Mean) diameter for this swarm
    real :: diameter

    ! Swarm mean relaxation time
    real :: tau

    ! Coefficient of restitution (1.0 - elastic, 0.0 - sticky)
    real :: rst

    ! Time step for the swarm
    real :: dt

    ! Eddy time interval for stochastic eddy interaction (SEIM) model
    integer :: time_eim

    ! Number of sub-steps; time sub-steps
    integer :: n_sub_steps

    ! Reflected and deposited particles on the walls and the escaped particles
    ! (These variables are declared real because it is currently not ...
    ! ... possible to store integer bnd-cell based variables in backup)
    real, allocatable :: n_reflected(:)
    real, allocatable :: n_deposited(:)
    real, allocatable :: n_escaped(:)
    integer           :: n_trapped

    ! Particle statistics
    logical :: statistics

    ! Logical array if cell has particles
    logical, allocatable :: cell_has_particles(:)

    ! Ensemble-averaged statistics for the swarm
    real,    allocatable :: u_mean(:), v_mean(:), w_mean(:)
    real,    allocatable :: uu(:), vv(:), ww(:), uv(:), uw(:), vw(:)
    integer, allocatable :: n_states(:)

    ! Gradients of flow modeled  quantity "zeta"
    ! (for SGS of Hybrid_Les_Rans model)
    real, allocatable :: v2_mod(:), v2_mod_x(:), v2_mod_y(:), v2_mod_z(:)

    ! SGS Brownian diffusion force
    real, allocatable :: f_fuka_x(:), f_fuka_y(:), f_fuka_z(:)

    ! Variable holding the subgrid scale (SGS)  model type
    ! (Must be part of type definition for multiple materials)
    integer :: subgrid_scale_model

    integer, allocatable :: i_work(:)
    logical, allocatable :: l_work(:)
    real,    allocatable :: r_work(:)

    ! Working arrays, buffers for parallel version
    ! (Keyword "parameter: not allowed inside a type
    ! declaration. One might think of making a function)
    integer :: N_I_VARS =  6
    integer :: N_L_VARS =  3
    integer :: N_R_VARS = 15

    contains
      procedure          :: Allocate_Swarm
      procedure, private :: Bounce_Particle
      procedure, private :: Move_Particle
      procedure, private :: Move_Trapped
      procedure, private :: Particle_Forces
      procedure, private :: Trap_Particle

  end type

  !---------------------!
  !   Model constants   !
  !---------------------!

  ! For Fukagata model
  real :: c_o   = 2.1
  real :: c_eps = 1.0

  contains

  ! Member procedures sorted alphabetically
  include 'Swarm_Mod/Advance_Particles.f90'
  include 'Swarm_Mod/Allocate_Swarm.f90'
  include 'Swarm_Mod/Bounce_Particle.f90'
  include 'Swarm_Mod/Calculate_Mean.f90'
  include 'Swarm_Mod/Check_Periodicity.f90'
  include 'Swarm_Mod/Exchange_Particles.f90'
  include 'Swarm_Mod/Grad_Modeled_Flow.f90'
  include 'Swarm_Mod/Move_Particle.f90'
  include 'Swarm_Mod/Move_Trapped.f90'
  include 'Swarm_Mod/Particle_Forces.f90'
  include 'Swarm_Mod/Particle_Time_Scale.f90'
  include 'Swarm_Mod/Print_Statistics.f90'
  include 'Swarm_Mod/Sgs_Discrete_Random_Walk.f90'
  include 'Swarm_Mod/Sgs_Fukagata.f90'
  include 'Swarm_Mod/Trap_Particle.f90'

  end module
