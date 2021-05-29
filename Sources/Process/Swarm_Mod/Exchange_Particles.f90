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
  type(Particle_Type), pointer :: Part
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
    swarm % l_work(:) = .false.
    swarm % r_work(:) = 0.0

    do k = 1, swarm % n_particles

      ! Take aliases for the particle
      Part => swarm % particle(k)

      !-----------------------------------------------------!
      !   Pack data for sending (all processors which ...   !
      !   ... send will put data in this globall pool)      !
      !-----------------------------------------------------!
      if(Part % proc .eq. this_proc) then
        i = (k-1) * swarm % N_I_VARS
        swarm % i_work(i + 1) = Part % proc  ! where it resides
        swarm % i_work(i + 2) = Part % buff  ! where it wants to go
        swarm % i_work(i + 3) = Part % cell
        swarm % i_work(i + 4) = Part % node
        swarm % i_work(i + 5) = Grid % comm % cell_glo(Part % cell)
        swarm % i_work(i + 6) = Grid % comm % node_glo(Part % node)

        i = (k-1) * swarm % N_L_VARS
        swarm % l_work(i + 1) = Part % deposited
        swarm % l_work(i + 2) = Part % escaped
        swarm % l_work(i + 3) = Part % trapped

        i = (k-1) * swarm % N_R_VARS
        swarm % r_work(i +  1) = Part % x_n
        swarm % r_work(i +  2) = Part % y_n
        swarm % r_work(i +  3) = Part % z_n
        swarm % r_work(i +  4) = Part % x_o
        swarm % r_work(i +  5) = Part % y_o
        swarm % r_work(i +  6) = Part % z_o
        swarm % r_work(i +  7) = Part % u
        swarm % r_work(i +  8) = Part % v
        swarm % r_work(i +  9) = Part % w
        swarm % r_work(i + 10) = Part % d
        swarm % r_work(i + 11) = Part % cfl
        swarm % r_work(i + 12) = Part % smooth_n
        swarm % r_work(i + 13) = Part % smooth_o
        swarm % r_work(i + 14) = Part % density
        swarm % r_work(i + 15) = Part % dens_fluid
      end if

    end do    ! through particles

    !-----------------------!
    !   Exchange the data   !
    !-----------------------!
    call Comm_Mod_Global_Sum_Int_Array (                   &
                  swarm % n_particles * swarm % N_I_VARS,  &
                  swarm % i_work)
    call Comm_Mod_Global_Lor_Log_Array (                   &
                  swarm % n_particles * swarm % N_L_VARS,  &
                  swarm % l_work)
    call Comm_Mod_Global_Sum_Real_Array(                   &
                  swarm % n_particles * swarm % N_R_VARS,  &
                  swarm % r_work)

    !-----------------------------------------!
    !   Distribute global data on particles   !
    !-----------------------------------------!
    do k = 1, swarm % n_particles

      ! Take alias
      Part => swarm % particle(k)

      i = (k-1) * swarm % N_I_VARS
      Part % proc = swarm % i_work(i + 1)
      Part % buff = swarm % i_work(i + 2)
      Part % cell = swarm % i_work(i + 3)
      Part % node = swarm % i_work(i + 4)

      ! Particle was not in this processor but wants to enter here
      if(Part % proc .ne. this_proc .and.  &
         Part % buff .eq. this_proc) then

        ! Set particle processor to correct value
        Part % proc = Part % buff

        ! Find the closest cell ...
        do c = 1, Grid % n_cells
          if(Grid % comm % cell_glo(c) .eq. swarm % i_work(i + 5)) then
            Part % cell = c
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
            Part % node = n
            goto 2
          end if
        end do
        print *, '# PANIC: Particle ', k, ' lost in transfer!'
        print *, '# Particle''s global node number is ', swarm % i_work(i + 6)
2       continue

      ! Particle was not in this processor and doesn't even want to ...
      ! ... enter here or particle was in this processor but has left it
      else if(Part % buff .ne. this_proc) then
        Part % proc = Part % buff
      end if

      i = (k-1) * swarm % N_L_VARS
      Part % deposited = swarm % l_work(i + 1)
      Part % escaped   = swarm % l_work(i + 2)
      Part % trapped   = swarm % l_work(i + 3)

      i = (k-1) * swarm % N_R_VARS
      Part % x_n        = swarm % r_work(i +  1)
      Part % y_n        = swarm % r_work(i +  2)
      Part % z_n        = swarm % r_work(i +  3)
      Part % x_o        = swarm % r_work(i +  4)
      Part % y_o        = swarm % r_work(i +  5)
      Part % z_o        = swarm % r_work(i +  6)
      Part % u          = swarm % r_work(i +  7)
      Part % v          = swarm % r_work(i +  8)
      Part % w          = swarm % r_work(i +  9)
      Part % d          = swarm % r_work(i + 10)
      Part % cfl        = swarm % r_work(i + 11)
      Part % smooth_n   = swarm % r_work(i + 12)
      Part % smooth_o   = swarm % r_work(i + 13)
      Part % density    = swarm % r_work(i + 14)
      Part % dens_fluid = swarm % r_work(i + 15)
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
