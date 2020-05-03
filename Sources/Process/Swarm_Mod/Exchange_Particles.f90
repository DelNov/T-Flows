!==============================================================================!
  subroutine Swarm_Mod_Exchange_Particles(swarm)
!------------------------------------------------------------------------------!
!   Advances all particles in the swarm.                                       !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Swarm_Type), target :: swarm
  integer                  :: k      ! particle number
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
        swarm % i_work(i + 3) = grid % comm % cell_glo(part % cell)

        i = (k-1) * N_R_VARS
        swarm % r_work(i +  1) = part % x_n
        swarm % r_work(i +  2) = part % y_n
        swarm % r_work(i +  3) = part % z_n
        swarm % r_work(i +  4) = part % u
        swarm % r_work(i +  5) = part % v
        swarm % r_work(i +  6) = part % w
        swarm % r_work(i +  7) = part % d
        swarm % r_work(i +  8) = part % cfl
        ! The following data is not really needed, at least not yet:
        ! swarm % r_work(i +  9) = part % density
        ! swarm % r_work(i + 10) = part % x_o
        ! swarm % r_work(i + 11) = part % y_o
        ! swarm % r_work(i + 12) = part % z_o
        ! swarm % r_work(i + 13) = part % rel_u
        ! swarm % r_work(i + 14) = part % rel_v
        ! swarm % r_work(i + 15) = part % rel_w
        ! swarm % r_work(i + 16) = part % rel_vel
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
      part % cell = swarm % i_work(i + 3)  ! holds global number for the moment

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
      part % x_n     = swarm % r_work(i +  1)
      part % y_n     = swarm % r_work(i +  2)
      part % z_n     = swarm % r_work(i +  3)
      part % u       = swarm % r_work(i +  4)
      part % v       = swarm % r_work(i +  5)
      part % w       = swarm % r_work(i +  6)
      part % d       = swarm % r_work(i +  7)
      part % cfl     = swarm % r_work(i +  8)
      ! The following data is not really needed, at least not yet:
      ! part % density = swarm % r_work(i +  9)
      ! part % x_o     = swarm % r_work(i + 10)
      ! part % y_o     = swarm % r_work(i + 11)
      ! part % z_o     = swarm % r_work(i + 12)
      ! part % rel_u   = swarm % r_work(i + 13)
      ! part % rel_v   = swarm % r_work(i + 14)
      ! part % rel_w   = swarm % r_work(i + 15)
      ! part % rel_vel = swarm % r_work(i + 16)
    end do

  end if

  end subroutine
