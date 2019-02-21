!==============================================================================!
  module Swarm_Mod
!------------------------------------------------------------------------------!
!   Module for Lagrangian particle tracking                                    !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod,  only: Grid_Type
  use Field_Mod, only: Field_Type
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
    real :: rho

    ! Particle surface area
    real :: area

    ! Particle Volume
    real :: vol

    ! Particle's diameter
    real :: d

    ! Particle time step
    real :: dtp

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
    type(Particle_Type), allocatable :: particles(:)

    ! Counter for depositing particles
    integer :: c_d

    ! Counter for escaping particles
    integer :: c_e

    ! Counter for reflected particles
    integer :: c_r

  end type

  contains

  include 'Swarm_Mod/Allocate_Particles.f90'
  include 'Swarm_Mod/Find_Nearest_Cell.f90'
  include 'Swarm_Mod/Find_Nearest_Node.f90'
  include 'Swarm_Mod/Find_Neighboring_Cells.f90'
  include 'Swarm_Mod/Interpolate_Velocity.f90'
  include 'Swarm_Mod/Particle_Forces.f90'

  end module
