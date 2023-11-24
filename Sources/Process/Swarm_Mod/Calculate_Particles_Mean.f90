!==============================================================================!
  subroutine Calculate_Particles_Mean(Swarm, k, n_stat_p)
!------------------------------------------------------------------------------!
!> This subroutine calculates time-averaged velocities and resolved Reynolds
!> stresses for particles within a swarm in T-Flows. It is integral for
!> analyzing particle dynamics in scale-resolving simulations such as LES.
!------------------------------------------------------------------------------!
! Functionality                                                                !
!                                                                              !
! * Time-averaged velocities: Computes mean velocities (u, v, w) for the       !
!   particles in each computational cell, crucial for ensemble averaging.      !
! * Resolved reynolds stresses: Calculates resolved Reynolds stresses (uu,     !
!   vv, ww, uv, uw, vw) for the particles, providing insights into turbulent   !
! * Tracking particle states: Maintains a count of the number of states        !
!   (time steps) for averaging, ensuring accurate statistical representation.  !
! * Compatibility with various models: Supports LES, DES, and other            !
!   scale-resolving turbulence models, adapting to the simulation context.     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Swarm_Type), target :: Swarm     !! the swarm of particles
  integer,       intent(in) :: k         !! particle number (rank)
  integer,       intent(in) :: n_stat_p  !! starting step for swarm statistics
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),     pointer :: Grid
  type(Field_Type),    pointer :: Flow
  type(Turb_Type),     pointer :: Turb
  type(Particle_Type), pointer :: Part
  integer                      :: c, l
  real                         :: m
!==============================================================================!

  if(.not. Swarm % statistics) return

  ! Take aliases
  Grid => Swarm % pnt_grid
  Flow => Swarm % pnt_flow
  Turb => Swarm % pnt_turb

  l = Time % Curr_Dt() - n_stat_p
  if(l > -1) then

    !---------------------------------!
    !   Scale-resolving simulations   !
    !---------------------------------!
    if(Turb % model .eq. LES_SMAGORINSKY    .or.  &
       Turb % model .eq. LES_DYNAMIC        .or.  &
       Turb % model .eq. LES_WALE           .or.  &
       Turb % model .eq. LES_TVM            .or.  &
       Turb % model .eq. HYBRID_LES_PRANDTL .or.  &
       Turb % model .eq. HYBRID_LES_RANS    .or.  &
       Turb % model .eq. DES_SPALART        .or.  &
       Turb % model .eq. DNS) then

      ! Take alias for the particle
      Part => Swarm % Particle(k)

      ! Cell in which the current particle resides
      c = Swarm % Particle(k) % cell

      ! Current number of states (for swarm quantity averaging) 
      m = real(Swarm % n_states(c))

      ! Mean velocities for swarm
      Swarm % u_mean(c) = (Swarm % u_mean(c) * m + Part % u) / (m+1)
      Swarm % v_mean(c) = (Swarm % v_mean(c) * m + Part % v) / (m+1)
      Swarm % w_mean(c) = (Swarm % w_mean(c) * m + Part % w) / (m+1)

      ! Resolved Reynolds stresses
      Swarm % uu(c) = (Swarm % uu(c) * m + Part % u * Part % u) / (m+1)
      Swarm % vv(c) = (Swarm % vv(c) * m + Part % v * Part % v) / (m+1)
      Swarm % ww(c) = (Swarm % ww(c) * m + Part % w * Part % w) / (m+1)

      Swarm % uv(c) = (Swarm % uv(c) * m + Part % u * Part % v) / (m+1)
      Swarm % uw(c) = (Swarm % uw(c) * m + Part % u * Part % w) / (m+1)
      Swarm % vw(c) = (Swarm % vw(c) * m + Part % v * Part % w) / (m+1)

      ! Increase the number of states of the cell
      Swarm % n_states(c) = Swarm % n_states(c) + 1

    end if
  end if

  end subroutine
