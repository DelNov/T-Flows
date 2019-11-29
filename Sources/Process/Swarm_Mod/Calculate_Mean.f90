!==============================================================================!
  subroutine Swarm_Mod_Calculate_Mean(swarm, k, n0, n1)
!------------------------------------------------------------------------------!
!   Calculates particle time averaged velocity                                 !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Turb_Type),  target :: turb
  type(Swarm_Type), target :: swarm
  integer                  :: n0, n1
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type), pointer :: flow
  type(Grid_Type),  pointer :: grid
  type(Particle_Type), pointer :: part
  integer                   :: c, n, sc, k
!==============================================================================!

  if(.not. swarm % swarm_statistics) return

  ! Take aliases
  grid => swarm % pnt_grid
  flow => swarm % pnt_flow

  if(n1 > n0) then

      !---------------------------------!
      !   Scale-resolving simulations   ! 
      !---------------------------------!
      if(turbulence_model .eq. LES_SMAGORINSKY    .or.  &
         turbulence_model .eq. LES_DYNAMIC        .or.  &
         turbulence_model .eq. LES_WALE           .or.  &
         turbulence_model .eq. HYBRID_LES_PRANDTL .or.  &
         turbulence_model .eq. DES_SPALART        .or.  &
         turbulence_model .eq. DNS) then

        ! Take alias for the particle
        part => swarm % particle(k)

        ! Cell in which the current particle resides
        c = swarm % particle(k) % cell

        n = swarm % n_states(c)

        ! Mean velocities
        swarm % u_mean(c) = (swarm % u_mean(c) * (1.*n) + part % u) / (1.*(n+1))
        swarm % v_mean(c) = (swarm % v_mean(c) * (1.*n) + part % v) / (1.*(n+1))
        swarm % w_mean(c) = (swarm % w_mean(c) * (1.*n) + part % w) / (1.*(n+1))

        ! Resolved Reynolds stresses
        swarm % uu(c) = (swarm % uu(c)*(1.*n) + part % u * part % u) / (1.*(n+1))
        swarm % vv(c) = (swarm % vv(c)*(1.*n) + part % v * part % v) / (1.*(n+1))
        swarm % ww(c) = (swarm % ww(c)*(1.*n) + part % w * part % w) / (1.*(n+1))

        swarm % uv(c) = (swarm % uv(c)*(1.*n) + part % u * part % v) / (1.*(n+1))
        swarm % uw(c) = (swarm % uw(c)*(1.*n) + part % u * part % w) / (1.*(n+1))
        swarm % vw(c) = (swarm % vw(c)*(1.*n) + part % v * part % w) / (1.*(n+1))

        ! Increaset the number of states for the cell
        swarm % n_states(c) = swarm % n_states(c) + 1
      end if
  end if

  end subroutine
