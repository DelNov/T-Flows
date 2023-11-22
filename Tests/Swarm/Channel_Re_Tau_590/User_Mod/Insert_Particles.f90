!==============================================================================!
  subroutine User_Mod_Insert_Particles(Flow, Turb, Vof, Swarm)
!------------------------------------------------------------------------------!
!   This function is called at the beginning of time step.                     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type), target     :: Flow
  type(Turb_Type),  target     :: Turb
  type(Vof_Type),   target     :: Vof
  type(Swarm_Type), target     :: Swarm
!----------------------------------[Locals]------------------------------------!
  type(Grid_Type), pointer :: Grid
  integer                  :: k, n_parts_in_buffers
  real                     :: c1, c2, c3   ! random variables
!------------------------------[Local parameters]------------------------------!
  real, parameter :: L1 = 6.28  ! streamwise
  real, parameter :: L2 = 3.14  ! spanwise
  real, parameter :: L3 = 2.0   ! wall-normal
!==============================================================================!

  ! Take alias(es)
  Grid => Flow % pnt_grid

  !-----------------------!
  !   24001st time step   !
  !-----------------------!
  if(Time % Curr_Dt() .eq. 24001) then  ! should be after the Flow is developed

    ! Track maximum number of particles
    Swarm % n_particles = Swarm % max_particles

    ! Browsing through all introduced particles
    do k = 1, Swarm % n_particles

      ! Generating random locations for particle
      call random_number(c1)
      call random_number(c2)
      call random_number(c3)
      call Swarm % Particle(k) % Insert_At(L1*c1, L2*c2, L3*c3,  &
                                           n_parts_in_buffers,   &
                                           Flow=Flow)
    end do
  end if

  end subroutine
