!==============================================================================!
  module Swarm_Mod
!------------------------------------------------------------------------------!
!   Module for Lagrangian particle tracking                                    !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod, only: HUGE, ONE_SIXTH, EARTH_G, PI
  use Grid_Mod,  only: Grid_Type
  use Var_Mod,   only: Var_Type
  use Field_Mod, only: Field_Type, density, viscosity
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !-------------------!
  !   Particle type   !
  !-------------------!
  type Particle_Type

    ! Particle's coordinates
    real :: x
    real :: y
    real :: z

    ! Particle's velocity
    real :: u
    real :: v
    real :: w

    ! Particle's density
    real :: density

    ! Particle surface area
    real :: area

    ! Particle Volume
    real :: vol

    ! Particle's diameter
    real :: d

    ! Particle time step
    real :: dt

    ! Particle relaxation time
    real :: tau

    ! The closest cell and node
    integer :: cell
    integer :: node

    ! Particle relative (magnitude) velocity
    real :: rel_vn  ! +ve value of relative velocity

    ! Particle relative velocity in y-dir (buoyancy force)
    real :: rel_vv
    real :: rel_vvn  ! +ve value of y-dir relative velocity

    ! Particle velocity magnitude
    real :: vel_p

    ! Particle Reynolds number (computed from relative velocity)
    real :: re

    ! Particle drag factor (from Re_p)
    real :: f    ! this is not to be confused with the drag coefficient

    ! Forces exerted on particle
    real :: fd                           ! drag force
    real :: fb                           ! buoyant force
    real :: ft                           ! total force

    ! Particle terminal velocity
    real :: vt

    ! Particle deposition and departure from domain 
    logical  :: deposited
    logical  :: escaped

    ! trapped BC type
    logical :: trapped

    ! Reflection BC type
    logical :: reflected

    ! Counter for testing
    integer :: counter

    ! Particle cfl number
    real :: cfl

    type(Grid_Type),  pointer :: pnt_grid  ! grid for which it is defined
    type(Field_Type), pointer :: pnt_flow  ! flow field for which it is defined

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

    ! Counter for depositing (d), escaped (e) and reflected (r) particles
    integer :: cnt_d
    integer :: cnt_e
    integer :: cnt_r

  end type

  contains

  include 'Swarm_Mod/Allocate_Particles.f90'
  include 'Swarm_Mod/Find_Nearest_Cell.f90'
  include 'Swarm_Mod/Find_Nearest_Node.f90'
  include 'Swarm_Mod/Find_Neighboring_Cells.f90'
  include 'Swarm_Mod/Interpolate_Velocity.f90'
  include 'Swarm_Mod/Particle_Forces.f90'

  end module
