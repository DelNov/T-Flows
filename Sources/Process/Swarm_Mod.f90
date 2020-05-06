!==============================================================================!
  module Swarm_Mod
!------------------------------------------------------------------------------!
!   Module for Lagrangian particle tracking                                    !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
  use Comm_Mod
  use Grid_Mod
  use Var_Mod
  use Field_Mod
  use Turb_Mod
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

    ! Particle Reynolds number (computed from relative velocity)
    real :: re

    ! Particle Courant number
    real :: cfl

    ! Particle drag factor (from Re_p)
    real :: f    ! this is not to be confused with the drag coefficient

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

  !----------------!
  !   Swarm type   !
  !----------------!
  type Swarm_Type

    type(Grid_Type),  pointer :: pnt_grid  ! grid for which it is defined
    type(Field_Type), pointer :: pnt_flow  ! flow field for which it is defined

    integer                          :: n_particles
    type(Particle_Type), allocatable :: particle(:)

    ! Density of this swarm
    real :: density

    ! (Mean) diameter for this swarm
    real :: diameter

    ! Coefficient of restitution (1.0 - elastic, 0.0 - sticky)
    real :: rst

    ! Time step for the swarm
    real :: dt

    ! Number of sub-steps; time sub-steps
    integer :: n_sub_steps

    ! Counter for depositing (d), escaped (e) and reflected (r) particles
    integer :: cnt_d
    integer :: cnt_e
    integer :: cnt_r

    ! Logical array if cell has particles
    logical, allocatable :: cell_has_particles(:)

    ! Working arrays, buffers for parallel version
    integer, allocatable :: i_work(:)
    logical, allocatable :: l_work(:)
    real,    allocatable :: r_work(:)

  end type

  integer, parameter   :: N_I_VARS = 3
  integer, parameter   :: N_L_VARS = 2
  integer, parameter   :: N_R_VARS = 8

  contains

  include 'Swarm_Mod/Advance_Particles.f90'
  include 'Swarm_Mod/Allocate.f90'
  include 'Swarm_Mod/Bounce_Particle.f90'
  include 'Swarm_Mod/Check_Periodicity.f90'
  include 'Swarm_Mod/Exchange_Particles.f90'
  include 'Swarm_Mod/Find_Nearest_Cell.f90'
  include 'Swarm_Mod/Find_Nearest_Node.f90'
  include 'Swarm_Mod/Move_Particle.f90'
  include 'Swarm_Mod/Particle_Forces.f90'

  end module
