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
  integer :: k, kk
  real    :: dx
!==============================================================================!

  !-------------------!
  !   1st time step   !
  !-------------------!
  if (n .eq. 1201) then     ! should be after the flow is developed

    ! Initializing both deposition and departure counters
    swarm % cnt_d = 0
    swarm % cnt_e = 0
    swarm % cnt_r = 0

    ! Allocate memory for each particle
    call Swarm_Mod_Allocate(flow, swarm)

    do k = 1, swarm % n_particles

      ! Placing particles (only at the 1st time step)
      kk = k-1
      dx = 20 * kk

      ! Take the diameter and density from the swarm
      swarm % particle(k) % d       = swarm % diameter
      swarm % particle(k) % density = swarm % density

      swarm % particle(k) % x = -0.00375 + dx * swarm % particle(k) % d
      swarm % particle(k) % y = 0.0599999
      swarm % particle(k) % z = 0.0

      ! Initializing particle's velocity
      swarm % particle(k) % u = 0
      swarm % particle(k) % v = 0
      swarm % particle(k) % w = 0

      ! Particle is in the domain (1st time step)
      swarm % particle(k) % deposited = .false.
      swarm % particle(k) % escaped   = .false.

      ! Define the type of BC at the wall for particles 
      ! (very important to set up)
      swarm % particle(k) % trapped   = .false.
      swarm % particle(k) % reflected = .true.

      ! Searching for the closest cell to start the algorithm (done once)
      call Swarm_Mod_Find_Nearest_Cell(swarm, k)

    end do

    print *, ""
    print *, "trapped particles =",swarm % cnt_d,  &
        ",","escaped particles =",swarm % cnt_e,       &
        ",","reflected particles =",swarm % cnt_r
  end if

  !----------------------!
  !   2nd time step on   !
  !----------------------!
  if (n .gt. 1201) then     ! should be after the flow is developed
    do k=1, swarm % n_particles

      call Swarm_Mod_Find_Neighboring_cells(swarm, k)

    end do

    print *, ""
    print *, "trapped particles =",swarm % cnt_d,  &
         ",","escaped particles =",swarm % cnt_e,       &
         ",","reflected particles =",swarm % cnt_r

  end if

  end subroutine
