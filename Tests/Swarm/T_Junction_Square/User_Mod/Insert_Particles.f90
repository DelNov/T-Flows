!==============================================================================!
  subroutine User_Mod_Insert_Particles(Flow, Turb, Vof, Swarm, n, time)
!------------------------------------------------------------------------------!
!   This function is called at the end of time step.                           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type), target :: Flow
  type(Turb_Type),  target :: Turb
  type(Vof_Type),   target :: Vof
  type(Swarm_Type), target :: Swarm
  integer, intent(in)      :: n         ! time step
  real,    intent(in)      :: time      ! physical time
!----------------------------------[Locals]------------------------------------!
  integer :: k, n_parts_in_buffers
  real    :: dx
!==============================================================================!

  !-------------------!
  !   1st time step   !
  !-------------------!
  if(n .eq. 1201) then     ! should be after the flow is developed

    ! Place the particles where you want them
    do k = 1, Swarm % n_particles

      ! Placing particles (only at the 1st time step)
      dx = 20 * (k - 1)

      Swarm % particle(k) % x_n = -0.00375 + dx * 2.5e-5
      Swarm % particle(k) % y_n = 0.0599999
      Swarm % particle(k) % z_n = 0.0

      Swarm % particle(k) % x_o = Swarm % particle(k) % x_n
      Swarm % particle(k) % y_o = Swarm % particle(k) % y_n
      Swarm % particle(k) % z_o = Swarm % particle(k) % z_n

      ! Searching for the closest cell and node to place the moved particle
      call Swarm % Particle(k) % Find_Nearest_Cell(n_parts_in_buffers)
      call Swarm % Particle(k) % Find_Nearest_Node()
    end do

  end if

  end subroutine
