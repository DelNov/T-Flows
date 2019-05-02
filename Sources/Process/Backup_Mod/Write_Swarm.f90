!==============================================================================!
  subroutine Backup_Mod_Write_Swarm(fh, disp, vc, swr)
!------------------------------------------------------------------------------!
!   Saves backup files name.backup                                             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer                  :: fh, disp, vc
  type(Swarm_Type), target :: swr
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),     pointer :: grid
  type(Particle_Type), pointer :: part
  integer                      :: i, k
!==============================================================================!

  ! Take aliases
  grid => swr % pnt_grid

  !--------------------------!
  !                          !
  !   Swarm (of particles)   !
  !                          !
  !--------------------------!

  ! Number of particles
  call Backup_Mod_Write_Int(fh, disp, vc, 'n_particles', swr % n_particles)

  ! Pack particle data in arrays
  do k = 1, swr % n_particles

    ! Take aliases for the particle
    part => swr % particle(k)

    i = (k-1) * n_i_vars
    i_work(i + 1) = part % proc  ! where it resides
    i_work(i + 2) = part % buff  ! where it wants to go
    i_work(i + 3) = grid % comm % cell_glo(part % cell)

    i = (k-1) * n_l_vars
    l_work(i + 1) = part % deposited
    l_work(i + 2) = part % escaped

    i = (k-1) * n_r_vars
    r_work(i + 1) = part % x_n
    r_work(i + 2) = part % y_n
    r_work(i + 3) = part % z_n
    r_work(i + 4) = part % u
    r_work(i + 5) = part % v
    r_work(i + 6) = part % w
    r_work(i + 7) = part % d
    r_work(i + 8) = part % cfl
  end do

  if(swr % n_particles > 0) then
    call Backup_Mod_Write_Int_Array(fh, disp, vc,         &
                                    'particle_int_data',  &
                                    i_work(1:n_i_vars*swr % n_particles))
    call Backup_Mod_Write_Log_Array(fh, disp, vc,         &
                                    'particle_log_data',  &
                                    l_work(1:n_l_vars*swr % n_particles))
    call Backup_Mod_Write_Real_Array(fh, disp, vc,          &
                                     'particle_real_data',  &
                                     r_work(1:n_r_vars*swr % n_particles))
  end if

  end subroutine
