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
  type(Comm_Type),     pointer :: Comm
  type(Particle_Type), pointer :: part
  integer                      :: i, k, n_part, n_parts_in_buffers
!==============================================================================!

  ! Take aliases
  Grid => swr  % pnt_grid
  Comm => Grid % Comm

  !--------------------------!
  !                          !
  !   Swarm (of particles)   !
  !                          !
  !--------------------------!

  ! Number of particles
  call Backup_Mod_Read_Int(Comm, fh, disp, vc, 'n_particles', n_part)

  swr % i_work(:) = 0
  swr % l_work(:) = .false.
  swr % r_work(:) = 0.0

  if(n_part > 0) then
    swr % n_particles = n_part
    call Backup_Mod_Read_Int_Array(Comm, fh, disp, vc,   &
                   'particle_int_data',                  &
                    swr % i_work(1 : swr % N_I_VARS*swr % n_particles))
    call Backup_Mod_Read_Log_Array(Comm, fh, disp, vc,   &
                   'particle_log_data',                  &
                    swr % l_work(1 : swr % N_L_VARS*swr % n_particles))
    call Backup_Mod_Read_Real_Array(Comm, fh, disp, vc,  &
                   'particle_real_data',                 &
                    swr % r_work(1 : swr % N_R_VARS*swr % n_particles))

    ! Pack Particle data in arrays
    do k = 1, swr % n_particles

      ! Take aliases for the Particle
      part => swr % Particle(k)

      i = (k-1) * swr % N_I_VARS
      part % proc = swr % i_work(i + 1)
      part % buff = swr % i_work(i + 2)
      part % cell = swr % i_work(i + 3)  ! holds global number for the moment

      i = (k-1) * swr % N_L_VARS
      part % deposited = swr % l_work(i + 1)
      part % escaped   = swr % l_work(i + 2)

      i = (k-1) * swr % N_R_VARS
      part % x_n  = swr % r_work(i + 1)
      part % y_n  = swr % r_work(i + 2)
      part % z_n  = swr % r_work(i + 3)
      part % u    = swr % r_work(i + 4)
      part % v    = swr % r_work(i + 5)
      part % w    = swr % r_work(i + 6)
      part % d    = swr % r_work(i + 7)
      part % cfl  = swr % r_work(i + 8)

    end do
  end if

  n_parts_in_buffers = 0
  do k = 1, swr % n_particles
    swr % Particle(k) % cell = 0
    swr % Particle(k) % node = 0
    swr % Particle(k) % proc = 0
    swr % Particle(k) % buff = 0
    call swr % Particle(k) % Find_Nearest_Cell(n_parts_in_buffers)
    call swr % Particle(k) % Find_Nearest_Node()
  end do

  end subroutine
