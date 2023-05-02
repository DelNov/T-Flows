!==============================================================================!
  subroutine User_Mod_Insert_Particles(Flow, Turb, Vof, Swarm)
!------------------------------------------------------------------------------!
!   This function is called at the end of time step.                           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type), target :: Flow
  type(Turb_Type),  target :: Turb
  type(Vof_Type),   target :: Vof
  type(Swarm_Type), target :: Swarm
!----------------------------------[Locals]------------------------------------!
  type(Grid_Type), pointer :: Grid
  integer                  :: k, n_parts_in_buffers, c
  real                     :: dx, rx, ry, rz
!==============================================================================!

  ! Take alias(es)
  Grid => Flow % pnt_grid

  !----------------------!
  !   1201st time step   !
  !----------------------!
  if(Time % Curr_Dt() .eq. 1201) then   ! should be after the flow is developed

    ! Compute all particles
    Swarm % n_particles = Swarm % max_particles

    ! Place the particles where you want them
    do k = 1, Swarm % n_particles
      dx = 20 * (k - 1)
      call Swarm % Particle(k) % Insert_At(-0.00375 + dx * 2.5e-5,  &
                                            0.0399999,              &
                                            0.0,                    &
                                            n_parts_in_buffers,     &
                                            Flow=Flow)
    end do

  end if

  end subroutine
