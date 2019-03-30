!==============================================================================!
  subroutine Backup_Mod_Read_Swarm(fh, disp, vc, swr)
!------------------------------------------------------------------------------!
!   Loads backup files name.backup                                             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer                  :: fh, disp, vc
  type(Swarm_Type), target :: swr
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),     pointer :: grid
  type(Particle_Type), pointer :: part
  integer                      :: i, k, c, n_part
!==============================================================================!

  ! Take aliases
  grid => swr % pnt_grid

  !--------------------------!
  !                          !
  !   Swarm (of particles)   !
  !                          !
  !--------------------------!

  ! Number of particles
  call Backup_Mod_Read_Int(fh, disp, vc, 'n_particles', n_part)

  if(n_part > 0) then
    swr % n_particles = n_part
    call Backup_Mod_Read_Int_Array(fh, disp, vc,         &
                                   'particle_int_data',  &
                                   i_work(1:n_i_vars*swr % n_particles))
    call Backup_Mod_Read_Log_Array(fh, disp, vc,         &
                                   'particle_log_data',  &
                                   l_work(1:n_l_vars*swr % n_particles))
    call Backup_Mod_Read_Real_Array(fh, disp, vc,          &
                                    'particle_real_data',  &
                                    r_work(1:n_r_vars*swr % n_particles))

    ! Pack particle data in arrays
    do k = 1, swr % n_particles

      ! Take aliases for the particle
      part => swr % particle(k)

      i = (k-1) * n_i_vars
      part % proc = i_work(i + 1)
      part % buff = i_work(i + 2)
      part % cell = i_work(i + 3)  ! holds global number for the moment

      i = (k-1) * n_l_vars
      part % deposited = l_work(i + 1)
      part % escaped   = l_work(i + 2)

      i = (k-1) * n_r_vars
      part % x_n  = r_work(i + 1)
      part % y_n  = r_work(i + 2)
      part % z_n  = r_work(i + 3)
      part % u    = r_work(i + 4)
      part % v    = r_work(i + 5)
      part % w    = r_work(i + 6)
      part % d    = r_work(i + 7)
      part % cfl  = r_work(i + 8)

      ! Searching for the closest cell and node to place the moved particle
      part % node = 0  ! force it to look for all cells
      call Swarm_Mod_Find_Nearest_Cell(swr, k)
      call Swarm_Mod_Find_Nearest_Node(swr, k)
    end do
  end if

  end subroutine
