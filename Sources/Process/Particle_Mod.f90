!==============================================================================!
  module Particle_Mod
!------------------------------------------------------------------------------!
!   Particle module used in Lagrangian particle tracking                       !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Vert_Mod
  use Field_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !-------------------!
  !   Particle type   !
  !-------------------!
  type, extends(Vert_Type) :: Particle_Type

    ! Flow field for which the particle is defined
    type(Field_Type), pointer :: pnt_flow

    ! Particle's velocity
    real :: u
    real :: v
    real :: w

    ! Particle's density
    real :: density

    ! Density of surrounding fluid
    real :: dens_fluid

    ! Particle's diameter
    real :: d

    ! Particle relaxation time
    real :: tau

    ! Particle relative velocity components and magnitude
    real :: rel_u
    real :: rel_v
    real :: rel_w
    real :: rel_vel

    ! Relative velocities (modeled flow quantity and particle's)
    real :: rel_u_mod
    real :: rel_v_mod
    real :: rel_w_mod

    ! Velocity fluctuations from stochastic eddy interaction (DRW model)
    real :: u_drw
    real :: v_drw
    real :: w_drw

    ! Particle Courant, Stokes and Reynolds numbers
    real :: cfl
    real :: st
    real :: re

    ! Particle drag factor (from Re_p)
    real :: f    ! this is not to be confused with the drag coefficient

    ! Particle terminal speed
    real :: vel_t

    ! Forces exerted on the particle
    real :: fd_x, fd_y, fd_z  ! drag force
    real :: fb_x, fb_y, fb_z  ! buoyant force

    ! New and old value of smoothed vof
    real :: smooth_n
    real :: smooth_o

    contains
      procedure :: Initialize_Particle
      procedure :: Insert_At

  end type

  contains
  include 'Particle_Mod/Initialize_Particle.f90'
  include 'Particle_Mod/Insert_At.f90'

  end module
