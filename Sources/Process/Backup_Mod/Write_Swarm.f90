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

  i_work(:) = 0
  l_work(:) = 0
  r_work(:) = 0.0

  ! Pack particle data in arrays
  do k = 1, swr % n_particles

    ! Take aliases for the particle
    part => swr % particle(k)

    if(part % proc .eq. this_proc) then
      i = (k-1) * N_I_VARS
      i_work(i + 1) = part % proc  ! where it resides
      i_work(i + 2) = part % buff  ! where it wants to go
      i_work(i + 3) = grid % comm % cell_glo(part % cell)

      i = (k-1) * N_L_VARS
      l_work(i + 1) = part % deposited
      l_work(i + 2) = part % escaped

      i = (k-1) * N_R_VARS
      r_work(i + 1) = part % x_n
      r_work(i + 2) = part % y_n
      r_work(i + 3) = part % z_n
      r_work(i + 4) = part % u
      r_work(i + 5) = part % v
      r_work(i + 6) = part % w
      r_work(i + 7) = part % d
      r_work(i + 8) = part % cfl
    end if  ! particle on this processor
  end do

  !-----------------------!
  !   Exchange the data   !
  !-----------------------!
  call Comm_Mod_Global_Sum_Int_Array (swr % n_particles * N_I_VARS, i_work)
  call Comm_Mod_Global_Sum_Real_Array(swr % n_particles * N_R_VARS, r_work)

  if(swr % n_particles > 0) then
    call Backup_Mod_Write_Int_Array(fh, disp, vc,         &
                                    'particle_int_data',  &
                                    i_work(1:N_I_VARS*swr % n_particles))
    call Backup_Mod_Write_Int_Array(fh, disp, vc,         &
                                    'particle_log_data',  &
                                    l_work(1:N_L_VARS*swr % n_particles))
    call Backup_Mod_Write_Real_Array(fh, disp, vc,          &
                                     'particle_real_data',  &
                                     r_work(1:N_R_VARS*swr % n_particles))
  end if

  end subroutine
