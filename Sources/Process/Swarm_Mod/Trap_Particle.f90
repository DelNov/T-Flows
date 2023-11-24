!==============================================================================!
  subroutine Trap_Particle(Swarm, k)
!------------------------------------------------------------------------------!
!> This subroutine is responsible for trapping particles at interfaces between
!> two phases within the computational domain. It is a key component in
!> simulations involving three-phase flows, two fluid phases resolved with VOF
!> and a solid phase, represented by particles.
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * Interface detection: Determines if a particle has crossed an interface   !
!     a critical aspect in three-phase flow simulations.                       !
!   * Particle trapping: Implements the logic to trap particles at interfaces, !
!     altering their state and behavior in response to the encountered surface.!
!   * Smoothed field evaluation: Utilizes the smoothed (VOF) field to          !
!     accurately assess the particle's position relative to the interface.     !
!   * Trapping flag management: Updates particle status flags to reflect       !
!     trapping, important for tracking particle dynamics and interactions.     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Swarm_Type), target :: Swarm  !! the swarm of particles
  integer, intent(in)       :: k      !! particle number (rank)
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),     pointer :: Grid
  type(Vof_Type),      pointer :: Vof
  type(Particle_Type), pointer :: Part
  logical,             pointer :: trapped              ! part. trapping flag
  integer                      :: c                    ! nearest cells, face
  real                         :: rx_n, ry_n, rz_n
!==============================================================================!

  ! Take aliases
  Grid      => Swarm % pnt_grid
  Vof       => Swarm % pnt_vof
  Part      => Swarm % Particle(k)
  trapped   => Part  % trapped

  c = Part % cell      ! index of the closest cell for interpolation

  ! Vector connecting new particle position and cell center
  rx_n = Part % x_n - Grid % xc(c)
  ry_n = Part % y_n - Grid % yc(c)
  rz_n = Part % z_n - Grid % zc(c)

  ! Value of smoothed vof at the new position
  Part % smooth_n = Vof % smooth % n(c)         &
                  + Vof % smooth % x(c) * rx_n  &
                  + Vof % smooth % y(c) * ry_n  &
                  + Vof % smooth % z(c) * rz_n

! PRINT *, Vof % Smooth % x(c), Vof % Smooth % y(c), Vof % Smooth % z(c)
! PRINT *, smooth_old, smooth_new

  ! Particle crossed an interface
  if( (Part % smooth_n - 0.5) * (Part % smooth_o - 0.5) < 0.0) then

    trapped = .true.
    Swarm % n_trapped = Swarm % n_trapped + 1

    print '(a,i6,a,1pe11.3,1pe11.3,1pe11.3)',    &
          ' # Particle ', k, ' trapped at  : ',  &
          Part % x_n, Part % y_n, Part % z_n

  end if  ! crossed an interface cell

  end subroutine
