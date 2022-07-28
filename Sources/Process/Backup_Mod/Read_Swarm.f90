!==============================================================================!
  subroutine Backup_Mod_Read_Swarm(disp, vc, Swr)
!------------------------------------------------------------------------------!
!   Loads backup files name.backup                                             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer(DP)              :: disp
  integer                  :: vc
  type(Swarm_Type), target :: Swr
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),     pointer :: Grid
  type(Comm_Type),     pointer :: Comm
  type(Particle_Type), pointer :: part
  integer                      :: i, k, n_part, n_parts_in_buffers
!==============================================================================!

  ! Take aliases
  Grid => Swr  % pnt_grid
  Comm => Grid % Comm

  !--------------------------!
  !                          !
  !   Swarm (of particles)   !
  !                          !
  !--------------------------!

  ! Number of particles
  call Backup_Mod_Read_Int(Comm, disp, vc, 'n_particles', n_part)

  Swr % i_work(:) = 0
  Swr % l_work(:) = .false.
  Swr % r_work(:) = 0.0

  if(n_part > 0) then
    Swr % n_particles = n_part
    call Backup_Mod_Read_Int_Array(Comm, disp, vc,   &
                   'particle_int_data',              &
                    Swr % i_work(1 : Swr % N_I_VARS*Swr % n_particles))
    call Backup_Mod_Read_Log_Array(Comm, disp, vc,   &
                   'particle_log_data',              &
                    Swr % l_work(1 : Swr % N_L_VARS*Swr % n_particles))
    call Backup_Mod_Read_Real_Array(Comm, disp, vc,  &
                   'particle_real_data',             &
                    Swr % r_work(1 : Swr % N_R_VARS*Swr % n_particles))

    ! Unpack particle data from arrays
    do k = 1, Swr % n_particles

      ! Take aliases for the Particle
      part => Swr % Particle(k)

      i = (k-1) * Swr % N_I_VARS
      part % proc = Swr % i_work(i + 1)
      part % buff = Swr % i_work(i + 2)
      part % cell = Swr % i_work(i + 3)  ! holds global number for the moment

      i = (k-1) * Swr % N_L_VARS
      part % deposited = Swr % l_work(i + 1)
      part % escaped   = Swr % l_work(i + 2)

      i = (k-1) * Swr % N_R_VARS
      part % x_n  = Swr % r_work(i + 1)
      part % y_n  = Swr % r_work(i + 2)
      part % z_n  = Swr % r_work(i + 3)
      part % u    = Swr % r_work(i + 4)
      part % v    = Swr % r_work(i + 5)
      part % w    = Swr % r_work(i + 6)
      part % d    = Swr % r_work(i + 7)
      part % cfl  = Swr % r_work(i + 8)

    end do
  end if

  n_parts_in_buffers = 0
  do k = 1, Swr % n_particles
    Swr % Particle(k) % cell = 0
    Swr % Particle(k) % node = 0
    Swr % Particle(k) % proc = 0
    Swr % Particle(k) % buff = 0
    call Swr % Particle(k) % Find_Nearest_Cell(n_parts_in_buffers)
    call Swr % Particle(k) % Find_Nearest_Node()
  end do

  end subroutine
