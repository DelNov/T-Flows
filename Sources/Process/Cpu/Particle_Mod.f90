!==============================================================================!
  module Particle_Mod
!------------------------------------------------------------------------------!
!> The Particle_Mod module is integral for Lagrangian particle tracking in
!> T-Flows, handled by module Swarm_Mod. Particle_Mod extends the capabilities
!> of the Vert_Type, adding specific attributes and methods tailored for
!> simulating and tracking particles within a flow field. This module is
!> crucial for analyzing particle dynamics in Swarm_Mod, including their
!> interaction with fluid forces and their movement in a computational domain.
!------------------------------------------------------------------------------!
! Features                                                                     !
!                                                                              !
! * Extension of Vert_Type: Inherits vertex properties and introduces          !
!   particle-specific attributes such as velocity, density, and diameter.      !
! * Flow field interaction: Links the particle with the flow field, allowing   !
!   for dynamic interaction and tracking within the flow field environment.    !
! * Velocity and forces: Manages particle's velocity components, drag, and     !
!   buoyant forces, providing detailed insights into particle motion.          !
! * Timescale and motion characteristics: Computes the particle's relaxation   !
!   time, Courant, Stokes, and Reynolds numbers, crucial for motion analysis.  !
! * Drag factor and terminal speed: Evaluates the drag factor based on         !
!   Reynolds number and determines particle's terminal speed.                  !
! * Stochastic interaction: Handles velocity fluctuations for particle-fluid   !
!   interaction models like the DRW model.                                     !
! * Initialization and insertion: Includes procedures to initialize and        !
!   position particles within the computational domain. !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Vof_Mod  ! contains Particle_Type
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !-------------------!
  !   Particle type   !
  !-------------------!
  !> Particle_Type extends the capabilities of Vert_Type for particle tracking,
  !> incorporating detailed attributes and methods specific to particle dynamics
  type, extends(Vert_Type) :: Particle_Type

    ! Flow field for which the particle is defined
    type(Field_Type), pointer :: pnt_flow  !! pointer to flow field

    ! Particle's velocity
    real :: u  !! particle's x-velocity component
    real :: v  !! particle's y-velocity component
    real :: w  !! particle's z-velocity component

    ! Particle's density
    real :: density  !! particle's density

    ! Density of surrounding fluid
    real :: dens_fluid  !! density of the surrounding fluid

    ! Particle's diameter
    real :: d  !! particle diameter

    ! Particle relaxation time
    real :: tau  !! particle relaxation time

    ! Particle relative velocity components and magnitude
    real :: rel_u    !! relative x-velocity component
    real :: rel_v    !! relative y-velocity component
    real :: rel_w    !! relative z-velocity component
    real :: rel_vel  !! magnitued of particle's relative velocity

    ! Relative velocities (modeled flow quantity and particle's)
    real :: rel_u_mod  !! relative x-velocity component (modeled)
    real :: rel_v_mod  !! relative y-velocity component (modeled)
    real :: rel_w_mod  !! relative z-velocity component (modeled)

    ! Velocity fluctuations from stochastic eddy interaction (DRW model)
    real :: u_drw  !! x-velocity fluctuation (DRW model)
    real :: v_drw  !! y-velocity fluctuation (DRW model)
    real :: w_drw  !! z-velocity fluctuation (DRW model)

    ! Particle Courant, Stokes and Reynolds numbers
    real :: cfl  !! particle's Courant number
    real :: st   !! particle's Stokes number
    real :: re   !! particle's Reynolds number

    ! Particle drag factor (from Re_p)
    real :: f    !! particle's drag factor (should not to
                 !! be confused with the drag coefficient)

    ! Particle terminal speed
    real :: vel_t  !! particle's terminal velocity

    ! Forces exerted on the particle
    real :: fd_x, fd_y, fd_z  !! drag force component
    real :: fb_x, fb_y, fb_z  !! buoyant force component

    ! New and old value of smoothed vof
    real :: smooth_n  !! new value of smoothed VOF
    real :: smooth_o  !! old value of smoothed VOF

    contains
      procedure :: Initialize_Particle
      procedure :: Insert_At

  end type

  contains

#   include "Particle_Mod/Initialize_Particle.f90"
#   include "Particle_Mod/Insert_At.f90"

  end module
