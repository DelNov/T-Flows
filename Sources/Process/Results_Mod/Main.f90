!==============================================================================!
  subroutine Results_Mod_Main(curr_dt, last_dt, time, n_dom,  &
                              flow, turb, mult, swarm, exit_now)
!------------------------------------------------------------------------------!
!   Main function for saving results (postprocessing and backup)               !
!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer               :: curr_dt, last_dt
  real                  :: time            ! physical time of the simulation
  integer               :: n_dom
  type(Field_Type)      :: flow(MD)        ! flow field
  type(Turb_Type)       :: turb(MD)        ! turbulence modelling
  type(Multiphase_Type) :: mult(MD)        ! multiphase modelling
  type(Swarm_Type)      :: swarm(MD)       ! swarm of particles
  logical               :: exit_now
!----------------------------------[Locals]------------------------------------!
  integer :: d         ! domain counter
  logical :: save_now
!==============================================================================!

  inquire(file='exit_now', exist=exit_now)
  inquire(file='save_now', exist=save_now)

  ! Is it time to save the backup file?
  if(curr_dt .eq. last_dt             .or.  &
     save_now                         .or.  &
     exit_now                         .or.  &
     Backup_Mod_Time_To_Save(curr_dt) .or.  &
     Info_Mod_Time_To_Exit()) then
    do d = 1, n_dom
      call Control_Mod_Switch_To_Domain(d)
      call Backup_Mod_Save(flow(d), swarm(d), turb(d), mult(d),  &
                           time, curr_dt, domain=d)
    end do
  end if

  ! Is it time to save results for post-processing?
  if(curr_dt .eq. last_dt              .or.  &
     save_now                          .or.  &
     exit_now                          .or.  &
     Results_Mod_Time_To_Save(curr_dt) .or.  &
     Info_Mod_Time_To_Exit()) then

    do d = 1, n_dom
      call Control_Mod_Switch_To_Domain(d)
      call Results_Mod_Save(flow(d), turb(d), mult(d), swarm(d), curr_dt,  &
                            plot_inside=.true., domain=d)
      call Results_Mod_Save(flow(d), turb(d), mult(d), swarm(d), curr_dt,  &
                            plot_inside=.false., domain=d)
      call Results_Mod_Save_Swarm(swarm(d), curr_dt)

      if(mult(d) % model .eq. VOLUME_OF_FLUID) then
      end if

      ! Write results in user-customized format
      call User_Mod_Save_Results(flow(d), turb(d), mult(d), swarm(d), curr_dt)
      call User_Mod_Save_Swarm(flow(d), turb(d), mult(d), swarm(d), curr_dt)

    end do  ! through domains
  end if

  if(save_now) then
    if(this_proc < 2) then
      open (9, file='save_now', status='old')
      close(9, status='delete')
    end if
  end if

  if(exit_now) then
    if(this_proc < 2) then
      open (9, file='exit_now', status='old')
      close(9, status='delete')
    end if
!   goto 2
  end if

  end subroutine
