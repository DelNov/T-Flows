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
  type(Grid_Type),     pointer :: Grid
  type(Comm_Type),     pointer :: Comm
  type(Particle_Type), pointer :: part
  integer                      :: i, k
!==============================================================================!

  ! Take aliases
  Grid => swr % pnt_grid
  Comm => Grid % Comm

  !--------------------------!
  !                          !
  !   Swarm (of particles)   !
  !                          !
  !--------------------------!

  ! Number of particles
  call Backup_Mod_Write_Int(Comm, fh, disp, vc, 'n_particles',  &
                                                swr % n_particles)

  swr % i_work(:) = 0
  swr % l_work(:) = .false.
  swr % r_work(:) = 0.0

  ! Pack particle data in arrays
  do k = 1, swr % n_particles

    ! Take aliases for the particle
    part => swr % particle(k)

    if(part % proc .eq. this_proc) then
      i = (k-1) * swr % N_I_VARS
      swr % i_work(i + 1) = part % proc  ! where it resides
      swr % i_work(i + 2) = part % buff  ! where it wants to go
      swr % i_work(i + 3) = Grid % comm % cell_glo(part % cell)

      i = (k-1) * swr % N_L_VARS
      swr % l_work(i + 1) = part % deposited
      swr % l_work(i + 2) = part % escaped

      i = (k-1) * swr % N_R_VARS
      swr % r_work(i + 1) = part % x_n
      swr % r_work(i + 2) = part % y_n
      swr % r_work(i + 3) = part % z_n
      swr % r_work(i + 4) = part % u
      swr % r_work(i + 5) = part % v
      swr % r_work(i + 6) = part % w
      swr % r_work(i + 7) = part % d
      swr % r_work(i + 8) = part % cfl
    end if  ! particle on this processor
  end do

  !-----------------------!
  !   Exchange the data   !
  !-----------------------!
  call Comm_Mod_Global_Sum_Int_Array (swr % n_particles * swr % N_I_VARS,  &
                                      swr % i_work)
  call Comm_Mod_Global_Sum_Real_Array(swr % n_particles * swr % N_R_VARS,  &
                                      swr % r_work)

  if(swr % n_particles > 0) then
    call Backup_Mod_Write_Int_Array(Comm, fh, disp, vc,   &
                   'particle_int_data',                   &
                    swr % i_work(1 : swr % N_I_VARS*swr % n_particles))
    call Backup_Mod_Write_Log_Array(Comm, fh, disp, vc,   &
                   'particle_log_data',                   &
                    swr % l_work(1 : swr % N_L_VARS*swr % n_particles))
    call Backup_Mod_Write_Real_Array(Comm, fh, disp, vc,  &
                   'particle_real_data',                  &
                    swr % r_work(1 : swr % N_R_VARS*swr % n_particles))
  end if

  end subroutine
