!==============================================================================!
  module Swarm_Mod
!------------------------------------------------------------------------------!
!   Module for Lagrangian particle tracking                                    !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod
  use Field_Mod
  use Turb_Mod
  use Control_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !-------------------!
  !   Particle type   !
  !-------------------!
  type Particle_Type

    ! Particle's coordinates; new and old
    real :: x_n, x_o
    real :: y_n, y_o
    real :: z_n, z_o

    ! Particle's velocity
    real :: u
    real :: v
    real :: w

    ! Particle's density
    real :: density

    ! Particle's diameter
    real :: d

    ! Particle relaxation time
    real :: tau

    ! The closest cell, node, boundary cell and face
    integer :: cell
    integer :: node
    integer :: bnd_cell
    integer :: bnd_face

    ! Particle relative velocity components and magnitude
    real :: rel_u
    real :: rel_v
    real :: rel_w
    real :: rel_vel

    ! Particle Reynolds number (computed from relative velocity and "flow-
    ! viscosity")
    real :: re

    ! Particle Courant number
    real :: cfl

    ! Particle drag factor (from Re_p)
    real :: f    ! this is not to be confused with the drag coefficient

    ! Particle terminal speed
    real :: vel_t

    ! Particle distance between two successive sub-time-steps
    real :: dsp 

    ! Forces exerted on the particle
    real :: fd_x, fd_y, fd_z  ! drag force
    real :: fb_x, fb_y, fb_z  ! buoyant force
    real :: ft_x, ft_y, ft_z  ! total force

    ! Particle deposition and departure from domain 
    logical :: deposited
    logical :: escaped

    ! Particle inside the subdomain
    integer :: proc
    integer :: buff

  end type

  ! Parameters describing turbulence model choice
  ! (Prime numbers starting from 30169  as a continuation for turb_mod)
  integer, parameter :: BROWNIAN_FUKAGATA  = 30169

  !----------------!
  !   Swarm type   !
  !----------------!
  type Swarm_Type

    type(Grid_Type),  pointer :: pnt_grid  ! grid for which it is defined
    type(Field_Type), pointer :: pnt_flow  ! flow field for which it is defined
    type(Turb_Type),  pointer :: pnt_turb  ! turb field for which it is defined

    integer                          :: n_particles
    type(Particle_Type), allocatable :: particle(:)

    ! Density of this swarm
    real :: density

    ! (Mean) diameter for this swarm
    real :: diameter

    ! Swarm mean relaxation time
    real :: tau

    ! Coefficient of restitution (1.0 - elastic, 0.0 - sticky)
    real :: rst

    ! Swarm's Stokes number
    real :: st

    ! Time step for the swarm
    real :: dt

    ! Number of sub-steps; time sub-steps
    integer :: n_sub_steps

    ! Counter for depositing (d), escaped (e) and reflected (r) particles
    integer :: cnt_d
    integer :: cnt_e
    integer :: cnt_r

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
    real, allocatable :: w_mod_x(:), w_mod_y(:), w_mod_z(:)

    ! Modeled TKE for swarm over the whole grid
    real, allocatable :: kin_mod(:)

    ! SGS Brownian diffusion force
    real, allocatable :: f_fuka_x(:), f_fuka_y(:), f_fuka_z(:)

    ! Counter (was initially added for sequential particle injection case)
    integer :: counter

    ! Variable holding the subgrid scale (SGS)  model type
    ! (Must be part of type definition for multiple materials)
    integer :: subgrid_scale_model

    integer, allocatable :: i_work(:)
    logical, allocatable :: l_work(:)
    real,    allocatable :: r_work(:)

  end type

  integer, parameter   :: N_I_VARS = 3
  integer, parameter   :: N_L_VARS = 2
  integer, parameter   :: N_R_VARS = 8

  !---------------------!
  !   Model constants   !
  !---------------------!

  ! For Fukagata model
  real :: c_o   = 2.1
  real :: c_eps = 1.0

  contains

  include 'Swarm_Mod/Advance_Particles.f90'
  include 'Swarm_Mod/Allocate.f90'
  include 'Swarm_Mod/Bounce_Particle.f90'
  include 'Swarm_Mod/Calculate_Mean.f90'
  include 'Swarm_Mod/Check_Periodicity.f90'
  include 'Swarm_Mod/Exchange_Particles.f90'
  include 'Swarm_Mod/Find_Nearest_Cell.f90'
  include 'Swarm_Mod/Find_Nearest_Node.f90'
  include 'Swarm_Mod/Grad_Modeled_Flow.f90'
  include 'Swarm_Mod/Move_Particle.f90'
  include 'Swarm_Mod/Particle_Forces.f90'
  include 'Swarm_Mod/Particle_Time_Scale.f90'
  include 'Swarm_Mod/Sgs_Fukagata.f90'

  end module
