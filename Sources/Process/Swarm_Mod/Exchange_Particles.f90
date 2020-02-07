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
  integer                      :: i, c
!==============================================================================!

  ! Take aliases for the swarm
  grid => swarm % pnt_grid

  !-----------------------------------------------!
  !   Exchange particle data between processors   !
  !-----------------------------------------------!
  if(n_proc > 1) then

    i_work(:) = 0
    r_work(:) = 0.0

    do k = 1, swarm % n_particles

      ! Take aliases for the particle
      part => swarm % particle(k)

      !-----------------------------------------------------!
      !   Pack data for sending (all processors which ...   !
      !   ... send will put data in this globall pool)      !
      !-----------------------------------------------------!
      if(part % proc .eq. this_proc) then
        i = (k-1) * N_I_VARS
        i_work(i + 1) = part % proc  ! where it resides
        i_work(i + 2) = part % buff  ! where it wants to go
        i_work(i + 3) = grid % comm % cell_glo(part % cell)

        i = (k-1) * N_R_VARS
        r_work(i +  1)  = part % x_n
        r_work(i +  2)  = part % y_n
        r_work(i +  3)  = part % z_n
        r_work(i +  4)  = part % u
        r_work(i +  5)  = part % v
        r_work(i +  6)  = part % w
        r_work(i +  7)  = part % d
        r_work(i +  8)  = part % cfl

        ! The following data is not really needed, at least not yet:
        ! r_work(i + 18) = part % density
        ! r_work(i + 19) = part % x_o
        ! r_work(i + 20) = part % y_o
        ! r_work(i + 21) = part % z_o
        ! r_work(i + 22) = part % rel_u
        ! r_work(i + 23) = part % rel_v
        ! r_work(i + 24) = part % rel_w
        ! r_work(i + 25) = part % rel_vel
      end if

    end do    ! through particles

    !-----------------------!
    !   Exchange the data   !
    !-----------------------!
    call Comm_Mod_Global_Sum_Int_Array (swarm % n_particles * N_I_VARS, i_work)
    call Comm_Mod_Global_Sum_Real_Array(swarm % n_particles * N_R_VARS, r_work)

    !-----------------------------------------!
    !   Distribute global data on particles   !
    !-----------------------------------------!
    do k = 1, swarm % n_particles

      ! Take alias
      part => swarm % particle(k)

      i = (k-1) * N_I_VARS
      part % proc = i_work(i + 1)
      part % buff = i_work(i + 2)
      part % cell = i_work(i + 3)  ! holds global number for the moment

      ! Particle was in this processor and wants to stay here
      if(part % proc .eq. this_proc .and. part % buff .eq. this_proc) then
        do c = 1, grid % n_cells
          if(grid % comm % cell_glo(c) .eq. part % cell) then
            part % cell = c
            exit
          end if
        end do
        call Swarm_Mod_Find_Nearest_Node(swarm, k)

      ! Particle was not in this processor but wants to enter here
      else if(part % proc .ne. this_proc .and. part % buff .eq. this_proc) then

        ! Set particle processor to correct value
        part % proc = part % buff

        ! If in its processor, ...
        if(part % proc .eq. this_proc) then

          ! ... find the closest cell ...
          do c = 1, grid % n_cells
            if(grid % comm % cell_glo(c) .eq. part % cell) then
              part % cell = c
              exit
            end if
          end do

          ! ... and the closest node.
          call Swarm_Mod_Find_Nearest_Node(swarm, k)
        else
          part % proc = 0
          part % buff = 0
        end if
      else
        part % proc = 0
        part % buff = 0
      end if

      i = (k-1) * N_R_VARS
      part % x_n     = r_work(i +  1)
      part % y_n     = r_work(i +  2)
      part % z_n     = r_work(i +  3)
      part % u       = r_work(i +  4)
      part % v       = r_work(i +  5)
      part % w       = r_work(i +  6)
      part % d       = r_work(i +  7)
      part % cfl     = r_work(i +  8)

      ! The following data is not really needed, at least not yet:
      ! part % density = r_work(i + 18)
      ! part % x_o     = r_work(i + 19)
      ! part % y_o     = r_work(i + 20)
      ! part % z_o     = r_work(i + 21)
      ! part % rel_u   = r_work(i + 22)
      ! part % rel_v   = r_work(i + 23)
      ! part % rel_w   = r_work(i + 24)
      ! part % rel_vel = r_work(i + 25)
    end do

    ! Refresh buffers for grid-base variables here
    ! (This is probably not needed, but it won't do harm)
    call Grid_Mod_Exchange_Real(grid, swarm % u_mean)
    call Grid_Mod_Exchange_Real(grid, swarm % v_mean)
    call Grid_Mod_Exchange_Real(grid, swarm % w_mean)
    call Grid_Mod_Exchange_Real(grid, swarm % uu)
    call Grid_Mod_Exchange_Real(grid, swarm % vv)
    call Grid_Mod_Exchange_Real(grid, swarm % ww)
    call Grid_Mod_Exchange_Real(grid, swarm % uv)
    call Grid_Mod_Exchange_Real(grid, swarm % uw)
    call Grid_Mod_Exchange_Real(grid, swarm % vw)

  end if

  end subroutine
