!==============================================================================!
  module Particle_Mod
!------------------------------------------------------------------------------!
!   Particle module used in Lagrangian particle tracking                       !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Point_Mod
  use Field_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !-------------------!
  !   Particle type   !
  !-------------------!
  type, extends(Point_Type) :: Particle_Type

    ! Flow field for which the particle is defined
    type(Field_Type), pointer :: pnt_flow

    ! Particle's old coordinates.  (New ones are in the parent.)
    real :: x_o
    real :: y_o
    real :: z_o

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

    ! Particle deposition and departure from domain
    ! (Should these be in the parent?  Time will tell.)
    logical :: deposited
    logical :: escaped

    contains
      procedure :: Initialize_Particle

  end type

  contains
  include 'Particle_Mod/Initialize_Particle.f90'

  end module
