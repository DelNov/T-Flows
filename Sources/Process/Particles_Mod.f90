!==============================================================================!
  module Particles_Mod
!------------------------------------------------------------------------------!
!   Module for Lagrangian particle tracking                                    !
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

    ! The closest cell
    integer :: cell

  end type

  contains

  include 'Particles_Mod/Find_Nearest_Cell.f90'
  include 'Particles_Mod/Interpolate_Velocity.f90'

  end module
