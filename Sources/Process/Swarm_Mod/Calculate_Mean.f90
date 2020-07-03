!==============================================================================!
  subroutine Swarm_Mod_Calculate_Mean(swarm, k, n, n_stat_p, ss)
!------------------------------------------------------------------------------!
!   Calculates particle time averaged velocity                                 !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Swarm_Type), target :: swarm
  integer                  :: k
  integer                  :: n
  integer                  :: n_stat_p
  integer                  :: ss
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),     pointer :: grid
  type(Field_Type),    pointer :: flow
  type(Turb_Type),     pointer :: turb
  type(Particle_Type), pointer :: part
  integer                      :: c, o, m, l
!==============================================================================!

  if(.not. swarm % statistics) return

  ! Take aliases
  grid => swarm % pnt_grid
  flow => swarm % pnt_flow
  turb => swarm % pnt_turb

  l = n - n_stat_p
  if(l > -1) then

    !---------------------------------!
    !   Scale-resolving simulations   !
    !---------------------------------!
    if(turb % model .eq. LES_SMAGORINSKY    .or.  &
       turb % model .eq. LES_DYNAMIC        .or.  &
       turb % model .eq. LES_WALE           .or.  &
       turb % model .eq. HYBRID_LES_PRANDTL .or.  &
       turb % model .eq. HYBRID_LES_RANS    .or.  &
       turb % model .eq. DES_SPALART        .or.  &
       turb % model .eq. DNS) then

      ! Take alias for the particle
      part => swarm % particle(k)

      ! Cell in which the current particle resides
      c = swarm % particle(k) % cell

      ! Current number of states (for swarm quantity averaging) 
      m = swarm % n_states(c)

      ! Current number of states (for particle velocity averaging) 
      n = part % n_states

      ! Mean velocities for "each" particle
      part % u_mean = (part % u_mean * (1.*o) + part % u) / (1.*(o+1))
      part % v_mean = (part % v_mean * (1.*o) + part % v) / (1.*(o+1))
      part % w_mean = (part % w_mean * (1.*o) + part % w) / (1.*(o+1))

      ! Mean velocities for swarm
      swarm % u_mean(c) = (swarm % u_mean(c) * (1.*m) + part % u) / (1.*(m+1))
      swarm % v_mean(c) = (swarm % v_mean(c) * (1.*m) + part % v) / (1.*(m+1))
      swarm % w_mean(c) = (swarm % w_mean(c) * (1.*m) + part % w) / (1.*(m+1))

      ! Resolved Reynolds stresses
      swarm % uu(c) = (swarm % uu(c)*(1.*m) + part % u * part % u) / (1.*(m+1))
      swarm % vv(c) = (swarm % vv(c)*(1.*m) + part % v * part % v) / (1.*(m+1))
      swarm % ww(c) = (swarm % ww(c)*(1.*m) + part % w * part % w) / (1.*(m+1))

      swarm % uv(c) = (swarm % uv(c)*(1.*m) + part % u * part % v) / (1.*(m+1))
      swarm % uw(c) = (swarm % uw(c)*(1.*m) + part % u * part % w) / (1.*(m+1))
      swarm % vw(c) = (swarm % vw(c)*(1.*m) + part % v * part % w) / (1.*(m+1))

      ! Increaset the number of states of the cell (for swarm) 
      swarm % n_states(c) = swarm % n_states(c) + 1

      ! Increaset the number of states for the particle
      part % n_states = part % n_states + 1

    end if
  end if

  end subroutine
