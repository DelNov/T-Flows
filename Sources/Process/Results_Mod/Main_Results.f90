!==============================================================================!
  subroutine Main_Results(Results, n_dom, Flow, Turb, Vof, Swarm, exit_now)
!------------------------------------------------------------------------------!
!   Main function for saving results (postprocessing and backup)               !
!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Results_Type) :: Results
  integer             :: n_dom
  type(Field_Type)    :: Flow(MD)        ! Flow field
  type(Turb_Type)     :: Turb(MD)        ! turbulence modelling
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
  if(Time % Curr_Dt() .eq. Time % Last_Dt() .or.  &
     save_now                               .or.  &
     exit_now                               .or.  &
     Backup % Time_To_Save_Backup()         .or.  &
     Info % Time_To_Exit()) then
    do d = 1, n_dom
      call Control % Switch_To_Domain(d)
      call Backup % Save(Flow(d), Turb(d), Vof(d), Swarm(d), dom=d)
    end do
  end if

  ! Is it time to save results for post-processing?
  if(Time % Curr_Dt() .eq. Time % Last_Dt() .or.  &
     save_now                               .or.  &
     exit_now                               .or.  &
     Results % Time_To_Save_Results()       .or.  &
     Info % Time_To_Exit()) then

    do d = 1, n_dom
      call Control % Switch_To_Domain(d)
      call Results % Save_Vtu_Results(Flow(d), Turb(d), Vof(d), Swarm(d),  &
                                      plot_inside=.true., domain=d)
      call Results % Save_Vtu_Results(Flow(d), Turb(d), Vof(d), Swarm(d),  &
                                      plot_inside=.false., domain=d)

      if(Flow(d) % with_interface) then
        if(Vof(d) % track_front) then
          call Results % Save_Vtu_Front(Vof(d) % Front)
        end if
        if(Vof(d) % track_surface) then
          call Results % Save_Vtu_Surf(Vof(d) % surf)
        end if
      end if

      ! Write results in user-customized format
      call User_Mod_Save_Results(Flow(d), Turb(d), Vof(d), Swarm(d), domain=d)
      if(Flow(d) % with_particles) then
        call User_Mod_Save_Swarm(Flow(d), Turb(d), Vof(d), Swarm(d), domain=d)
      end if

    end do  ! through domains
  end if

  ! Is it time to save particles for post-processing?
  if(Time % Curr_Dt() .eq. Time % Last_Dt() .or.  &
     save_now                               .or.  &
     exit_now                               .or.  &
     Results % Time_To_Save_Swarm()         .or.  &
     Info % Time_To_Exit()) then

    do d = 1, n_dom
      call Control % Switch_To_Domain(d)
      if(Flow(d) % with_particles) then
        call Results % Save_Vtu_Swarm(Swarm(d), domain=d)

        call User_Mod_Save_Swarm(Flow(d), Turb(d), Vof(d), Swarm(d), domain=d)
      end if

    end do  ! through domains
  end if

  if(save_now) then
    if(First_Proc()) then
      open (9, file='save_now', status='old')
      close(9, status='delete')
    end if
  end if

  if(exit_now) then
    if(First_Proc()) then
      open (9, file='exit_now', status='old')
      close(9, status='delete')
    end if
  end if

  end subroutine
