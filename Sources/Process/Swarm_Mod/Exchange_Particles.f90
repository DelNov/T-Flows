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
  type(Grid_Type),     pointer :: grid
  type(Particle_Type), pointer :: part
  integer                      :: i, c, n
!==============================================================================!

  ! Take aliases for the swarm
  grid => swarm % pnt_grid

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
        i = (k-1) * N_I_VARS
        swarm % i_work(i + 1) = part % proc  ! where it resides
        swarm % i_work(i + 2) = part % buff  ! where it wants to go
        swarm % i_work(i + 3) = part % cell
        swarm % i_work(i + 4) = part % node
        swarm % i_work(i + 5) = grid % comm % cell_glo(part % cell)
        swarm % i_work(i + 6) = grid % comm % node_glo(part % node)

        i = (k-1) * N_R_VARS
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
    call Comm_Mod_Global_Sum_Int_Array (swarm % n_particles * N_I_VARS,  &
                                        swarm % i_work)
    call Comm_Mod_Global_Sum_Real_Array(swarm % n_particles * N_R_VARS,  &
                                        swarm % r_work)

    !-----------------------------------------!
    !   Distribute global data on particles   !
    !-----------------------------------------!
    do k = 1, swarm % n_particles

      ! Take alias
      part => swarm % particle(k)

      i = (k-1) * N_I_VARS
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
        do c = 1, grid % n_cells
          if(grid % comm % cell_glo(c) .eq. swarm % i_work(i + 5)) then
            part % cell = c
            goto 1
          end if
        end do
        print *, '# PANIC: Particle ', k, ' lost in transfer!'
        print *, '# Particle''s global cell number is ', swarm % i_work(i + 5)
1       continue

        ! ... and the closest node.
        do n = 1, grid % n_nodes
          if(grid % comm % node_glo(n) .eq. swarm % i_work(i + 6)) then
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

      i = (k-1) * N_R_VARS

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

    !---------------------------------------------------------!
    !                                                         !
    !   Refresh buffers for grid-base variables here          !
    !   (This is probably not needed, but it won't do harm)   !
    !                                                         !
    !---------------------------------------------------------!
    if(swarm % statistics) then
      call Grid_Mod_Exchange_Cells_Int (grid, swarm % n_states)
      call Grid_Mod_Exchange_Cells_Real(grid, swarm % u_mean)
      call Grid_Mod_Exchange_Cells_Real(grid, swarm % v_mean)
      call Grid_Mod_Exchange_Cells_Real(grid, swarm % w_mean)
      call Grid_Mod_Exchange_Cells_Real(grid, swarm % uu)
      call Grid_Mod_Exchange_Cells_Real(grid, swarm % vv)
      call Grid_Mod_Exchange_Cells_Real(grid, swarm % ww)
      call Grid_Mod_Exchange_Cells_Real(grid, swarm % uv)
      call Grid_Mod_Exchange_Cells_Real(grid, swarm % uw)
      call Grid_Mod_Exchange_Cells_Real(grid, swarm % vw)
    end if

  end if

  end subroutine
