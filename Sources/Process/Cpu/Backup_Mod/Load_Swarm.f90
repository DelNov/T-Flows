!==============================================================================!
  subroutine Load_Swarm(Backup, Swarm, disp, vc)
!------------------------------------------------------------------------------!
!   Loads backup files name.backup                                             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Backup_Type)       :: Backup
  type(Swarm_Type), target :: Swarm
  integer(DP)              :: disp
  integer                  :: vc
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),     pointer :: Grid
  type(Comm_Type),     pointer :: Comm
  type(Particle_Type), pointer :: part
  integer                      :: i, k, n_part, n_parts_in_buffers
!==============================================================================!

  ! Take aliases
  Grid => Swarm  % pnt_grid
  Comm => Grid % Comm

  !--------------------------!
  !                          !
  !   Swarm (of particles)   !
  !                          !
  !--------------------------!

  ! Number of particles
  call Backup % Load_Int(Comm, disp, vc, 'n_particles', n_part)

  Swarm % i_work(:) = 0
  Swarm % l_work(:) = .false.
  Swarm % r_work(:) = 0.0

  if(n_part > 0) then
    Swarm % n_particles = n_part
    call Backup % Load_Int_Array(Comm, disp, vc,    &
                  'particle_int_data',              &
                  Swarm % i_work(1 : Swarm % N_I_VARS*Swarm % n_particles))
    call Backup % Load_Log_Array(Comm, disp, vc,    &
                  'particle_log_data',              &
                  Swarm % l_work(1 : Swarm % N_L_VARS*Swarm % n_particles))
    call Backup % Load_Real_Array(Comm, disp, vc,   &
                  'particle_real_data',             &
                  Swarm % r_work(1 : Swarm % N_R_VARS*Swarm % n_particles))

    ! Unpack particle data from arrays
    do k = 1, Swarm % n_particles

      ! Take aliases for the Particle
      part => Swarm % Particle(k)

      i = (k-1) * Swarm % N_I_VARS
      part % proc = Swarm % i_work(i + 1)
      part % buff = Swarm % i_work(i + 2)
      part % cell = Swarm % i_work(i + 3)  ! holds global number for the moment

      i = (k-1) * Swarm % N_L_VARS
      part % deposited = Swarm % l_work(i + 1)
      part % escaped   = Swarm % l_work(i + 2)

      i = (k-1) * Swarm % N_R_VARS
      part % x_n  = Swarm % r_work(i + 1)
      part % y_n  = Swarm % r_work(i + 2)
      part % z_n  = Swarm % r_work(i + 3)
      part % u    = Swarm % r_work(i + 4)
      part % v    = Swarm % r_work(i + 5)
      part % w    = Swarm % r_work(i + 6)
      part % d    = Swarm % r_work(i + 7)
      part % cfl  = Swarm % r_work(i + 8)

    end do
  end if

  n_parts_in_buffers = 0
  do k = 1, Swarm % n_particles
    Swarm % Particle(k) % cell = 0
    Swarm % Particle(k) % node = 0
    Swarm % Particle(k) % proc = 0
    Swarm % Particle(k) % buff = 0
    call Swarm % Particle(k) % Find_Nearest_Cell(n_parts_in_buffers)
    call Swarm % Particle(k) % Find_Nearest_Node()
  end do

  end subroutine
