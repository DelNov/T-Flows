#include "../../Shared/Assert.h90"
#include "../../Shared/Browse.h90"
#include "../../Shared/Macros.h90"
#include "../../Shared/Unused.h90"

!==============================================================================!
  module Swarm_Mod
!------------------------------------------------------------------------------!
!> The Swarm_Mod module is dedicated to Lagrangian particle tracking in
!> T-Flows. It manages swarms of particles, providing mechanisms for their
!> interaction with flow fields, turbulence modeling, and surface interactions.
!------------------------------------------------------------------------------!
! Features and Functionalities                                                 !
!                                                                              !
! * Particle integration: Manages the movement and interaction of particles    !
!   with the fluid and boundaries within a defined computational grid.         !
! * Turbulence interaction: Incorporates turbulence models like Brownian       !
!   motion (Fukagata) and Discrete Random Walk (DRW) for particle tracking.    !
! * Swarm dynamics: Handles the dynamic properties of particle swarms such     !
!   as density, diameter, relaxation time, and statistics.                     !
! * Particle-surface interactions: Manages particle behavior upon collision    !
!   with surfaces, including bouncing, trapping, and deposition.               !
! * Parallel processing: Supports particle exchange and tracking in parallel   !
!   computational environments, ensuring accurate swarm dynamics.              !
! * Statistical analysis: Provides functionality for calculating and           !
!   printing ensemble-averaged statistics of particle motion and interaction.  !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Particle_Mod
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
  !> Defines the properties and behaviors of a particle swarm within the
  !> computational grid. Extends particle dynamics to collective behavior
  !> and interaction with the flow field, turbulence and interfaces
  type Swarm_Type

    type(Grid_Type),  pointer :: pnt_grid  !! grid for which it is defined
    type(Field_Type), pointer :: pnt_flow  !! flow field for which it is defined
    type(Turb_Type),  pointer :: pnt_turb  !! turb field for which it is defined
    type(Vof_Type),   pointer :: pnt_vof   !! vof, for three-phase situations

    integer                          :: max_particles = 0  ! read from control
    integer                          :: n_particles   = 0
    type(Particle_Type), allocatable :: Particle(:)

    ! Density of this swarm
    real :: density  !! density of this swarm

    ! (Mean) diameter for this swarm
    real :: diameter  !! diameter in this swarm

    ! Swarm mean relaxation time
    real :: tau  !! mean relaxation time for the swarm

    ! Coefficient of restitution (1.0 - elastic, 0.0 - sticky)
    real :: rst  !! coefficient of restitution (1.0 - elastic, 0.0 - sticky)

    ! Time step for the swarm
    real :: dt  !! time step for the swarm

    ! Eddy time interval for stochastic eddy interaction (SEIM) model
    integer :: time_eim  !! eddy time interval for stochastic
                         !! eddy interaction (SEIM) model
    ! Number of sub-steps; time sub-steps
    integer :: n_sub_steps  !! number of sub-steps in time

    ! Reflected and deposited particles on the walls and the escaped particles
    ! (These variables are declared real because it is currently not ...
    ! ... possible to store integer bnd-cell based variables in backup)
    real, allocatable :: n_reflected(:)  !! number of reflected particles
    real, allocatable :: n_deposited(:)  !! number of deposited particles
    real, allocatable :: n_escaped(:)    !! number of escaped particles
    integer           :: n_trapped       !! number of trapped particles, in
                                         !! three-phase flow situations
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
    real, allocatable :: v2_mod(:)    !! modeled v2 quantity from zeta-f model
    real, allocatable :: v2_mod_x(:)  !! gradient of modeled v2 in x-direction
    real, allocatable :: v2_mod_y(:)  !! gradient of modeled v2 in x-direction
    real, allocatable :: v2_mod_z(:)  !! gradient of modeled v2 in x-direction

    ! SGS Brownian diffusion force
    real, allocatable :: f_fuka_x(:), f_fuka_y(:), f_fuka_z(:)

    ! Variable holding the subgrid scale (SGS)  model type
    ! (Must be part of type definition for multiple materials)
    integer :: subgrid_scale_model

    integer, allocatable :: i_work(:)  !! integer work array
    logical, allocatable :: l_work(:)  !! logical work array
    real,    allocatable :: r_work(:)  !! real work array

    ! Working arrays, buffers for parallel version
    ! (Keyword "parameter: not allowed inside a type
    ! declaration. One might think of making a function)
    integer :: N_I_VARS =  6  !! number of integer arrays
    integer :: N_L_VARS =  3  !! number of logical arrays
    integer :: N_R_VARS = 15  !! number of real arrays

    contains
      procedure          :: Advance_Particles
      procedure, private :: Bounce_Particle
      procedure, private :: Calculate_Particles_Mean
      procedure, private :: Check_Periodicity
      procedure          :: Create_Swarm
      procedure          :: Exchange_Particles
      procedure, private :: Grad_Modeled_Flow
      procedure, private :: Move_Particle
      procedure, private :: Move_Trapped
      procedure, private :: Particle_Forces
      procedure, private :: Particle_Time_Scale
      procedure, private :: Print_Swarm_Statistics
      procedure, private :: Sgs_Discrete_Random_Walk
      procedure, private :: Sgs_Fukagata
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
#   include "Swarm_Mod/Advance_Particles.f90"
#   include "Swarm_Mod/Create_Swarm.f90"
#   include "Swarm_Mod/Bounce_Particle.f90"
#   include "Swarm_Mod/Calculate_Particles_Mean.f90"
#   include "Swarm_Mod/Check_Periodicity.f90"
#   include "Swarm_Mod/Exchange_Particles.f90"
#   include "Swarm_Mod/Grad_Modeled_Flow.f90"
#   include "Swarm_Mod/Move_Particle.f90"
#   include "Swarm_Mod/Move_Trapped.f90"
#   include "Swarm_Mod/Particle_Forces.f90"
#   include "Swarm_Mod/Particle_Time_Scale.f90"
#   include "Swarm_Mod/Print_Swarm_Statistics.f90"
#   include "Swarm_Mod/Sgs_Discrete_Random_Walk.f90"
#   include "Swarm_Mod/Sgs_Fukagata.f90"
#   include "Swarm_Mod/Trap_Particle.f90"

  end module
