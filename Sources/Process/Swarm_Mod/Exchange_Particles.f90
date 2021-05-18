!==============================================================================!
  subroutine Swarm_Mod_Exchange_Particles(swarm)
!------------------------------------------------------------------------------!
!   Advances all particles in the swarm.                                       !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Swarm_Type), target     :: swarm
  integer                      :: k      ! particle number
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),     pointer :: Grid
  type(Particle_Type), pointer :: part
  integer                      :: i, c, n
!==============================================================================!

  ! Take aliases for the swarm
  Grid => swarm % pnt_grid

  !-----------------------------------------------!
  !                                               !
  !   Exchange particle data between processors   !
  !                                               !
  !-----------------------------------------------!
  if(n_proc > 1) then

    swarm % i_work(:) = 0
    swarm % r_work(:) = 0.0

    do k = 1, swarm % n_particles

      ! Take aliases for the particle
      part => swarm % particle(k)

      !-----------------------------------------------------!
      !   Pack data for sending (all processors which ...   !
      !   ... send will put data in this globall pool)      !
      !-----------------------------------------------------!
      if(part % proc .eq. this_proc) then
        i = (k-1) * swarm % N_I_VARS
        swarm % i_work(i + 1) = part % proc  ! where it resides
        swarm % i_work(i + 2) = part % buff  ! where it wants to go
        swarm % i_work(i + 3) = part % cell
        swarm % i_work(i + 4) = part % node
        swarm % i_work(i + 5) = Grid % comm % cell_glo(part % cell)
        swarm % i_work(i + 6) = Grid % comm % node_glo(part % node)

        i = (k-1) * swarm % N_R_VARS
        swarm % r_work(i +  1) = part % x_n
        swarm % r_work(i +  2) = part % y_n
        swarm % r_work(i +  3) = part % z_n
        swarm % r_work(i +  4) = part % x_o
        swarm % r_work(i +  5) = part % y_o
        swarm % r_work(i +  6) = part % z_o
        swarm % r_work(i +  7) = part % u
        swarm % r_work(i +  8) = part % v
        swarm % r_work(i +  9) = part % w
        swarm % r_work(i + 10) = part % d
        swarm % r_work(i + 11) = part % cfl
      end if

    end do    ! through particles

    !-----------------------!
    !   Exchange the data   !
    !-----------------------!
    call Comm_Mod_Global_Sum_Int_Array (                   &
                  swarm % n_particles * swarm % N_I_VARS,  &
                  swarm % i_work)
    call Comm_Mod_Global_Sum_Real_Array(                   &
                  swarm % n_particles * swarm % N_R_VARS,  &
                  swarm % r_work)

    !-----------------------------------------!
    !   Distribute global data on particles   !
    !-----------------------------------------!
    do k = 1, swarm % n_particles

      ! Take alias
      part => swarm % particle(k)

      i = (k-1) * swarm % N_I_VARS
      part % proc = swarm % i_work(i + 1)
      part % buff = swarm % i_work(i + 2)
      part % cell = swarm % i_work(i + 3)
      part % node = swarm % i_work(i + 4)

      ! Particle was not in this processor but wants to enter here
      if(part % proc .ne. this_proc .and.  &
         part % buff .eq. this_proc) then

        ! Set particle processor to correct value
        part % proc = part % buff

        ! Find the closest cell ...
        do c = 1, Grid % n_cells
          if(Grid % comm % cell_glo(c) .eq. swarm % i_work(i + 5)) then
            part % cell = c
            goto 1
          end if
        end do
        print *, '# PANIC: Particle ', k, ' lost in transfer!'
        print *, '# This might be because the mesh is too coarse or ' // &
                 'number of sub timesteps is too small'
        print *, '# Particle''s global cell number is ', swarm % i_work(i + 5)
1       continue

        ! ... and the closest node.
        do n = 1, Grid % n_nodes
          if(Grid % comm % node_glo(n) .eq. swarm % i_work(i + 6)) then
            part % node = n
            goto 2
          end if
        end do
        print *, '# PANIC: Particle ', k, ' lost in transfer!'
        print *, '# Particle''s global node number is ', swarm % i_work(i + 6)
2       continue

      ! Particle was not in this processor and doesn't even want to ...
      ! ... enter here or particle was in this processor but has left it
      else if(part % buff .ne. this_proc) then
        part % proc = part % buff
      end if

      i = (k-1) * swarm % N_R_VARS

      part % x_n = swarm % r_work(i +  1)
      part % y_n = swarm % r_work(i +  2)
      part % z_n = swarm % r_work(i +  3)
      part % x_o = swarm % r_work(i +  4)
      part % y_o = swarm % r_work(i +  5)
      part % z_o = swarm % r_work(i +  6)
      part % u   = swarm % r_work(i +  7)
      part % v   = swarm % r_work(i +  8)
      part % w   = swarm % r_work(i +  9)
      part % d   = swarm % r_work(i + 10)
      part % cfl = swarm % r_work(i + 11)
    end do

    !-------------------------------------------------------!
    !                                                       !
    !   Refresh buffers for Grid-base variables here        !
    !   (This is probably only needed for post-processing   !
    !    if buffers are plotted as well, but it is fine.)   !
    !                                                       !
    !-------------------------------------------------------!
    call Grid % Exchange_Cells_Real(swarm % n_reflected)
    call Grid % Exchange_Cells_Real(swarm % n_deposited)
    call Grid % Exchange_Cells_Real(swarm % n_escaped)

    if(swarm % statistics) then
      call Grid % Exchange_Cells_Int (swarm % n_states)
      call Grid % Exchange_Cells_Real(swarm % u_mean)
      call Grid % Exchange_Cells_Real(swarm % v_mean)
      call Grid % Exchange_Cells_Real(swarm % w_mean)
      call Grid % Exchange_Cells_Real(swarm % uu)
      call Grid % Exchange_Cells_Real(swarm % vv)
      call Grid % Exchange_Cells_Real(swarm % ww)
      call Grid % Exchange_Cells_Real(swarm % uv)
      call Grid % Exchange_Cells_Real(swarm % uw)
      call Grid % Exchange_Cells_Real(swarm % vw)
    end if

  end if

  end subroutine
