!==============================================================================!
  subroutine Backup_Mod_Write_Swarm(disp, vc, Swr)
!------------------------------------------------------------------------------!
!   Saves backup files name.backup                                             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer(DP)              :: disp
  integer                  :: vc
  type(Swarm_Type), target :: Swr
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),     pointer :: Grid
  type(Comm_Type),     pointer :: Comm
  type(Particle_Type), pointer :: Part
  integer                      :: i, k
!==============================================================================!

  ! Take aliases
  Grid => Swr % pnt_grid
  Comm => Grid % Comm

  !--------------------------!
  !                          !
  !   Swarm (of particles)   !
  !                          !
  !--------------------------!

  ! Number of particles
  call Backup_Mod_Write_Int(Comm, disp, vc, 'n_particles',  &
                                            Swr % n_particles)

  !----------------------------------------------!
  !   Write only if there are active particles   !
  !----------------------------------------------!
  if(Swr % n_particles > 0) then

    Swr % i_work(:) = 0
    Swr % l_work(:) = .false.
    Swr % r_work(:) = 0.0

    ! Pack particle data in arrays
    do k = 1, Swr % n_particles

      ! Take aliases for the particle
      Part => Swr % Particle(k)

      if(Part % proc .eq. this_proc) then
        i = (k-1) * Swr % N_I_VARS
        Swr % i_work(i + 1) = Part % proc  ! where it resides
        Swr % i_work(i + 2) = Part % buff  ! where it wants to go
        Swr % i_work(i + 3) = Grid % comm % cell_glo(Part % cell)

        i = (k-1) * Swr % N_L_VARS
        Swr % l_work(i + 1) = Part % deposited
        Swr % l_work(i + 2) = Part % escaped

        i = (k-1) * Swr % N_R_VARS
        Swr % r_work(i + 1) = Part % x_n
        Swr % r_work(i + 2) = Part % y_n
        Swr % r_work(i + 3) = Part % z_n
        Swr % r_work(i + 4) = Part % u
        Swr % r_work(i + 5) = Part % v
        Swr % r_work(i + 6) = Part % w
        Swr % r_work(i + 7) = Part % d
        Swr % r_work(i + 8) = Part % cfl
      end if  ! particle on this processor
    end do

    !-----------------------!
    !   Exchange the data   !
    !-----------------------!
    call Comm_Mod_Global_Sum_Int_Array (Swr % n_particles * Swr % N_I_VARS,  &
                                        Swr % i_work)
    call Comm_Mod_Global_Sum_Real_Array(Swr % n_particles * Swr % N_R_VARS,  &
                                        Swr % r_work)

    call Backup_Mod_Write_Int_Array(Comm, disp, vc,   &
                   'particle_int_data',               &
                    Swr % i_work(1 : Swr % N_I_VARS*Swr % n_particles))
    call Backup_Mod_Write_Log_Array(Comm, disp, vc,   &
                   'particle_log_data',               &
                    Swr % l_work(1 : Swr % N_L_VARS*Swr % n_particles))
    call Backup_Mod_Write_Real_Array(Comm, disp, vc,  &
                   'particle_real_data',              &
                    Swr % r_work(1 : Swr % N_R_VARS*Swr % n_particles))
  end if

  end subroutine
