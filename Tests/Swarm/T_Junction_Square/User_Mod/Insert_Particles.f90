!==============================================================================!
  subroutine User_Mod_Insert_Particles(flow, turb, Vof, swarm, n, time)
!------------------------------------------------------------------------------!
!   This function is called at the end of time step.                           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type), target :: flow
  type(Turb_Type),  target :: turb
  type(Vof_Type),   target :: Vof
  type(Swarm_Type), target :: swarm
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
    do k = 1, swarm % n_particles

      ! Placing particles (only at the 1st time step)
      dx = 20 * (k - 1)

      swarm % particle(k) % x_n = -0.00375 + dx * 2.5e-5
      swarm % particle(k) % y_n = 0.0599999
      swarm % particle(k) % z_n = 0.0

      swarm % particle(k) % x_o = swarm % particle(k) % x_n
      swarm % particle(k) % y_o = swarm % particle(k) % y_n
      swarm % particle(k) % z_o = swarm % particle(k) % z_n

      ! Searching for the closest cell and node to place the moved particle
      call Swarm_Mod_Find_Nearest_Cell(swarm, k, n_parts_in_buffers)
      call Swarm_Mod_Find_Nearest_Node(swarm, k)
    end do

  end if

  end subroutine
