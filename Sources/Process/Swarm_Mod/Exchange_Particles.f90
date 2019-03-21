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
  type(Grid_Type),     pointer     :: grid
  type(Particle_Type), pointer     :: part
  integer                          :: i, c
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

      ! Pack data for sending (all processors which ...
      ! ... send will put data in this globall pool)
      if(part % proc .eq. this_proc) then
        i = (k-1) * n_i_vars
        i_work(i + 1) = part % proc  ! where it resides
        i_work(i + 2) = part % buff  ! where it wants to go
        i_work(i + 3) = grid % comm % cell_glo(part % cell)

        i = (k-1) * n_r_vars
        r_work(i +  1) = part % x_n
        r_work(i +  2) = part % y_n
        r_work(i +  3) = part % z_n
        r_work(i +  4) = part % u
        r_work(i +  5) = part % v
        r_work(i +  6) = part % w
        r_work(i +  7) = part % d
        r_work(i +  8) = part % cfl
        ! The following data is not really needed, at least not yet:
        ! r_work(i +  9) = part % density
        ! r_work(i + 10) = part % x_o
        ! r_work(i + 11) = part % y_o
        ! r_work(i + 12) = part % z_o
        ! r_work(i + 13) = part % rel_u
        ! r_work(i + 14) = part % rel_v
        ! r_work(i + 15) = part % rel_w
        ! r_work(i + 16) = part % rel_vel
      end if

    end do    ! through particles

    ! Exchange the data
    call Comm_Mod_Global_Sum_Int_Array (i_work, swarm % n_particles * n_i_vars)
    call Comm_Mod_Global_Sum_Real_Array(r_work, swarm % n_particles * n_r_vars)

    ! Distribute global data on particles
    do k = 1, swarm % n_particles

      ! Take alias
      part => swarm % particle(k)

      i = (k-1) * n_i_vars
      part % proc = i_work(i + 1)
      part % buff = i_work(i + 2)
      part % cell = i_work(i + 3)  ! holds global number for the moment

      ! If particle crossed processor boundary
      if(part % buff .ne. part % proc) then

        ! Set particle processor to correct value
        part % proc = part % buff

        ! If in its processor, ...
        if(part % proc .eq. this_proc) then

          ! ... find the closest cell ...
          do c = 1, grid % n_cells
            if(grid % comm % cell_glo(c) .eq. part % cell) then
              part % cell = c
            end if
          end do

          ! ... and the closest node.
          call Swarm_Mod_Find_Nearest_Node(swarm, k)

        end if
      end if

      i = (k-1) * n_r_vars
      part % x_n     = r_work(i +  1)
      part % y_n     = r_work(i +  2)
      part % z_n     = r_work(i +  3)
      part % u       = r_work(i +  4)
      part % v       = r_work(i +  5)
      part % w       = r_work(i +  6)
      part % d       = r_work(i +  7)
      part % cfl     = r_work(i +  8)
      ! The following data is not really needed, at least not yet:
      ! part % density = r_work(i +  9)
      ! part % x_o     = r_work(i + 10)
      ! part % y_o     = r_work(i + 11)
      ! part % z_o     = r_work(i + 12)
      ! part % rel_u   = r_work(i + 13)
      ! part % rel_v   = r_work(i + 14)
      ! part % rel_w   = r_work(i + 15)
      ! part % rel_vel = r_work(i + 16)
    end do

  end if

  end subroutine
