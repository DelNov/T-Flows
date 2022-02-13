!==============================================================================!
  subroutine Main_Results(Results,                           &
                          curr_dt, last_dt, time, n_dom,     &
                          Flow, turb, Vof, Swarm, exit_now)
!------------------------------------------------------------------------------!
!   Main function for saving results (postprocessing and backup)               !
!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Results_Type) :: Results
  integer             :: curr_dt, last_dt
  real                :: time            ! physical time of the simulation
  integer             :: n_dom
  type(Field_Type)    :: Flow(MD)        ! Flow field
  type(Turb_Type)     :: turb(MD)        ! turbulence modelling
  type(Vof_Type)      :: Vof(MD)         ! multiphase modelling with vof
  type(Swarm_Type)    :: Swarm(MD)       ! swarm of particles
  logical             :: exit_now
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
      call Backup_Mod_Save(Flow(d), Swarm(d), turb(d), Vof(d),  &
                           time, curr_dt, domain=d)
    end do
  end if

  ! Is it time to save results for post-processing?
  if(curr_dt .eq. last_dt            .or.  &
     save_now                        .or.  &
     exit_now                        .or.  &
     Results % Time_To_Save(curr_dt) .or.  &
     Info_Mod_Time_To_Exit()) then

    do d = 1, n_dom
      call Control_Mod_Switch_To_Domain(d)
      call Results % Save_Results(Flow(d), turb(d), Vof(d), Swarm(d),  &
                                  curr_dt, plot_inside=.true., domain=d)
      call Results % Save_Results(Flow(d), turb(d), Vof(d), Swarm(d),  &
                                  curr_dt, plot_inside=.false., domain=d)
      call Results % Save_Swarm(Swarm(d), curr_dt)

      if(Flow(d) % with_interface) then
        if(Vof(d) % track_front) then
          call Results % Save_Front(Vof(d) % Front, curr_dt)
        end if
        if(Vof(d) % track_surface) then
          call Results % Save_Surf(Vof(d) % surf, curr_dt)
        end if
      end if

      ! Write results in user-customized format
      call User_Mod_Save_Results(Flow(d), turb(d), Vof(d), Swarm(d), curr_dt)
      call User_Mod_Save_Swarm(Flow(d), turb(d), Vof(d), Swarm(d), curr_dt)

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
  end if

  end subroutine
