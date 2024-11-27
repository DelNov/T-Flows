!==============================================================================!
  subroutine Save_Swarm(Backup, Swarm, disp, vc)
!------------------------------------------------------------------------------!
!   Saves backup files name.backup                                             !
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
  type(Particle_Type), pointer :: Part
  integer                      :: i, k
!==============================================================================!

  ! Take aliases
  Grid => Swarm % pnt_grid
  Comm => Grid % Comm

  !--------------------------!
  !                          !
  !   Swarm (of particles)   !
  !                          !
  !--------------------------!

  ! Number of particles
  call Backup % Save_Int(Comm, disp, vc, 'n_particles', Swarm % n_particles)

  !----------------------------------------------!
  !   Write only if there are active particles   !
  !----------------------------------------------!
  if(Swarm % n_particles > 0) then

    Swarm % i_work(:) = 0
    Swarm % l_work(:) = .false.
    Swarm % r_work(:) = 0.0

    ! Pack particle data in arrays
    do k = 1, Swarm % n_particles

      ! Take aliases for the particle
      Part => Swarm % Particle(k)

      if(Part % proc .eq. This_Proc()) then
        i = (k-1) * Swarm % N_I_VARS
        Swarm % i_work(i + 1) = Part % proc  ! where it resides
        Swarm % i_work(i + 2) = Part % buff  ! where it wants to go
        Swarm % i_work(i + 3) = Grid % Comm % cell_glo(Part % cell)

        i = (k-1) * Swarm % N_L_VARS
        Swarm % l_work(i + 1) = Part % deposited
        Swarm % l_work(i + 2) = Part % escaped

        i = (k-1) * Swarm % N_R_VARS
        Swarm % r_work(i + 1) = Part % x_n
        Swarm % r_work(i + 2) = Part % y_n
        Swarm % r_work(i + 3) = Part % z_n
        Swarm % r_work(i + 4) = Part % u
        Swarm % r_work(i + 5) = Part % v
        Swarm % r_work(i + 6) = Part % w
        Swarm % r_work(i + 7) = Part % d
        Swarm % r_work(i + 8) = Part % cfl
      end if  ! particle on this processor
    end do

    !-----------------------!
    !   Exchange the data   !
    !-----------------------!
    call Global % Sum_Int_Array (Swarm % n_particles*Swarm % N_I_VARS,  &
                                 Swarm % i_work)
    call Global % Sum_Real_Array(Swarm % n_particles*Swarm % N_R_VARS,  &
                                 Swarm % r_work)

    call Backup % Save_Int_Array(Comm, disp, vc,   &
                'particle_int_data',               &
                Swarm % i_work(1 : Swarm % N_I_VARS*Swarm % n_particles))
    call Backup % Save_Log_Array(Comm, disp, vc,   &
                'particle_log_data',               &
                Swarm % l_work(1 : Swarm % N_L_VARS*Swarm % n_particles))
    call Backup % Save_Real_Array(Comm, disp, vc,  &
                'particle_real_data',              &
                Swarm % r_work(1 : Swarm % N_R_VARS*Swarm % n_particles))
  end if

  end subroutine
