!==============================================================================!
  subroutine Exchange_Particles(Swarm)
!------------------------------------------------------------------------------!
!>  This subroutine handles the exchange of particle data between processors
!>  in a parallel computing environment. It is essential for ensuring the
!>  consistency and accuracy of particle tracking across processor boundaries.
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * Data preparation: Organizes particle data for exchange, including        !
!     processor and buffer information, position, velocity, and other          !
!     particle properties.                                                     !
!   * Data exchange: Facilitates the sharing of particle data among different  !
!     processors, ensuring that each processor has the correct information     !
!     about particles that move across boundaries.                             !
!   * Particle redistribution: Updates particles’ position and state based on  !
!     the received data, handling cases where particles enter or leave a       !
!     processor's sub-domain.                                                  !
!   * Grid-based variable update: Refreshes grid-based variables like particle !
!     statistics and flow quantities, which are necessary for accurate         !
!     particle tracking and post-processing.                                   !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Swarm_Type), target :: Swarm  !! the swarm of particles
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),     pointer :: Grid
  type(Particle_Type), pointer :: Part
  integer                      :: i, k, c, n
!==============================================================================!

  ! Take aliases for the Swarm
  Grid => Swarm % pnt_grid

  !-----------------------------------------------!
  !                                               !
  !   Exchange particle data between processors   !
  !                                               !
  !-----------------------------------------------!
  if(Parallel_Run()) then

    Swarm % i_work(:) = 0
    Swarm % l_work(:) = .false.
    Swarm % r_work(:) = 0.0

    do k = 1, Swarm % n_particles

      ! Take aliases for the particle
      Part => Swarm % Particle(k)

      !-----------------------------------------------------!
      !   Pack data for sending (all processors which ...   !
      !   ... send will put data in this globall pool)      !
      !-----------------------------------------------------!
      if(Part % proc .eq. This_Proc()) then
        i = (k-1) * Swarm % N_I_VARS
        Swarm % i_work(i + 1) = Part % proc  ! where it resides
        Swarm % i_work(i + 2) = Part % buff  ! where it wants to go
        Swarm % i_work(i + 3) = Part % cell
        Swarm % i_work(i + 4) = Part % node
        Swarm % i_work(i + 5) = Grid % Comm % cell_glo(Part % cell)
        Swarm % i_work(i + 6) = Grid % Comm % node_glo(Part % node)

        i = (k-1) * Swarm % N_L_VARS
        Swarm % l_work(i + 1) = Part % deposited
        Swarm % l_work(i + 2) = Part % escaped
        Swarm % l_work(i + 3) = Part % trapped

        i = (k-1) * Swarm % N_R_VARS
        Swarm % r_work(i +  1) = Part % x_n
        Swarm % r_work(i +  2) = Part % y_n
        Swarm % r_work(i +  3) = Part % z_n
        Swarm % r_work(i +  4) = Part % x_o
        Swarm % r_work(i +  5) = Part % y_o
        Swarm % r_work(i +  6) = Part % z_o
        Swarm % r_work(i +  7) = Part % u
        Swarm % r_work(i +  8) = Part % v
        Swarm % r_work(i +  9) = Part % w
        Swarm % r_work(i + 10) = Part % d
        Swarm % r_work(i + 11) = Part % cfl
        Swarm % r_work(i + 12) = Part % smooth_n
        Swarm % r_work(i + 13) = Part % smooth_o
        Swarm % r_work(i + 14) = Part % density
        Swarm % r_work(i + 15) = Part % dens_fluid
      end if

    end do    ! through particles

    !-----------------------!
    !   Exchange the data   !
    !-----------------------!
    call Global % Sum_Int_Array (                          &
                  Swarm % n_particles * Swarm % N_I_VARS,  &
                  Swarm % i_work)
    call Global % Lor_Log_Array (                          &
                  Swarm % n_particles * Swarm % N_L_VARS,  &
                  Swarm % l_work)
    call Global % Sum_Real_Array(                          &
                  Swarm % n_particles * Swarm % N_R_VARS,  &
                  Swarm % r_work)

    !-----------------------------------------!
    !   Distribute global data on particles   !
    !-----------------------------------------!
    do k = 1, Swarm % n_particles

      ! Take alias
      Part => Swarm % Particle(k)

      i = (k-1) * Swarm % N_I_VARS
      Part % proc = Swarm % i_work(i + 1)
      Part % buff = Swarm % i_work(i + 2)
      Part % cell = Swarm % i_work(i + 3)
      Part % node = Swarm % i_work(i + 4)

      ! Particle was not in this processor but wants to enter here
      if(Part % proc .ne. This_Proc() .and.  &
         Part % buff .eq. This_Proc()) then

        ! Set particle processor to correct value
        Part % proc = Part % buff

        ! Find the closest cell ...
        do c = 1, Grid % n_cells
          if(Grid % Comm % cell_glo(c) .eq. Swarm % i_work(i + 5)) then
            Part % cell = c
            goto 1
          end if
        end do
        print *, '# PANIC: Particle ', k, ' lost in transfer!'
        print *, '# This might be because the mesh is too coarse or ' // &
                 'number of sub timesteps is too small'
        print *, '# Particle''s global cell number is ', Swarm % i_work(i + 5)
1       continue

        ! ... and the closest node.
        do n = 1, Grid % n_nodes
          if(Grid % Comm % node_glo(n) .eq. Swarm % i_work(i + 6)) then
            Part % node = n
            goto 2
          end if
        end do
        print *, '# PANIC: Particle ', k, ' lost in transfer!'
        print *, '# Particle''s global node number is ', Swarm % i_work(i + 6)
2       continue

      ! Particle was not in this processor and doesn't even want to ...
      ! ... enter here or particle was in this processor but has left it
      else if(Part % buff .ne. This_Proc()) then
        Part % proc = Part % buff
      end if

      i = (k-1) * Swarm % N_L_VARS
      Part % deposited = Swarm % l_work(i + 1)
      Part % escaped   = Swarm % l_work(i + 2)
      Part % trapped   = Swarm % l_work(i + 3)

      i = (k-1) * Swarm % N_R_VARS
      Part % x_n        = Swarm % r_work(i +  1)
      Part % y_n        = Swarm % r_work(i +  2)
      Part % z_n        = Swarm % r_work(i +  3)
      Part % x_o        = Swarm % r_work(i +  4)
      Part % y_o        = Swarm % r_work(i +  5)
      Part % z_o        = Swarm % r_work(i +  6)
      Part % u          = Swarm % r_work(i +  7)
      Part % v          = Swarm % r_work(i +  8)
      Part % w          = Swarm % r_work(i +  9)
      Part % d          = Swarm % r_work(i + 10)
      Part % cfl        = Swarm % r_work(i + 11)
      Part % smooth_n   = Swarm % r_work(i + 12)
      Part % smooth_o   = Swarm % r_work(i + 13)
      Part % density    = Swarm % r_work(i + 14)
      Part % dens_fluid = Swarm % r_work(i + 15)
    end do

    !-------------------------------------------------------!
    !                                                       !
    !   Refresh buffers for Grid-base variables here        !
    !   (This is probably only needed for post-processing   !
    !    if buffers are plotted as well, but it is fine.)   !
    !                                                       !
    !-------------------------------------------------------!
    call Grid % Exchange_Cells_Real(Swarm % n_reflected)
    call Grid % Exchange_Cells_Real(Swarm % n_deposited)
    call Grid % Exchange_Cells_Real(Swarm % n_escaped)

    if(Swarm % statistics) then
      call Grid % Exchange_Cells_Int (Swarm % n_states)
      call Grid % Exchange_Cells_Real(Swarm % u_mean)
      call Grid % Exchange_Cells_Real(Swarm % v_mean)
      call Grid % Exchange_Cells_Real(Swarm % w_mean)
      call Grid % Exchange_Cells_Real(Swarm % uu)
      call Grid % Exchange_Cells_Real(Swarm % vv)
      call Grid % Exchange_Cells_Real(Swarm % ww)
      call Grid % Exchange_Cells_Real(Swarm % uv)
      call Grid % Exchange_Cells_Real(Swarm % uw)
      call Grid % Exchange_Cells_Real(Swarm % vw)
    end if

  end if

  end subroutine
