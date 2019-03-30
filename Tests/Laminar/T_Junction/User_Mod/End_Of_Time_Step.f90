!==============================================================================!
  subroutine User_Mod_End_Of_Time_Step(flow, swarm, n, time)
!------------------------------------------------------------------------------!
!   This function is called at the end of time step.                           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type), target :: flow
  type(Swarm_Type), target :: swarm
  integer                  :: n     ! time step
  real                     :: time  ! physical time
!----------------------------------[Locals]------------------------------------!
  integer :: k
  real    :: dx
!==============================================================================!

  !-------------------!
  !   1st time step   !
  !-------------------!
  if(n .eq. 1201) then     ! should be after the flow is developed

    ! Initializing both deposition and departure counters
    swarm % cnt_d = 0
    swarm % cnt_e = 0
    swarm % cnt_r = 0

    ! Place the particles where you want them
    do k = 1, swarm % n_particles

      ! Placing particles (only at the 1st time step)
      dx = 20 * (k - 1)

      swarm % particle(k) % x_n = -0.00375 + dx * 2.5e-5
      swarm % particle(k) % y_n = 0.0599999
      swarm % particle(k) % z_n = 0.0

      ! Searching for the closest cell and node to place the moved particle
      call Swarm_Mod_Find_Nearest_Cell(swarm, k)
      call Swarm_Mod_Find_Nearest_Node(swarm, k)
    end do

  end if

  !----------------------!
  !   2nd time step on   !
  !----------------------!
  if(n .gt. 1201) then     ! should be after the flow is developed
    call Swarm_Mod_Advance_Particles(swarm)
  end if

  if(this_proc < 2) then
    write(*,'(a,i4,a,i4,a,i4)'),                       &
             "# trapped particles: ",  swarm % cnt_d,  &
             " escaped particles: ",   swarm % cnt_e,  &
             " reflected particles: ", swarm % cnt_r
  end if

  end subroutine
