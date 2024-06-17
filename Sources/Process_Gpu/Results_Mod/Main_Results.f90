!==============================================================================!
  subroutine Main_Results(Results, Turb, Flow, Grid, n_dom)
!------------------------------------------------------------------------------!
!>  The Main_Results subroutine serves as the central function for saving
!>  results in T-Flows. It orchestrates the process of saving numerical results
!>  for post-processing and backups, ensuring that data is outputted correctly
!>  across different domains and under various conditions.
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * Conditional result saving: Determines whether to save initial conditions !
!     and handles exit conditions for result saving.                           !
!   * Domain-specific processing: Iterates through domains to save results     !
!     and backups for each domain separately.                                  !
!   * Post-processing output: Saves results in .vtu format for visualization,  !
!     including detailed data like flow fields, turbulence modeling data,      !
!     VOF multiphase modeling, and particle swarms.
!   * User-customized output: Allows for integration of user-defined result    !
!     saving routines through calls to User_Mod_Save_Fields and                !
!     User_Mod_Save_Swarm.                                                     !
!   * Handles file-based triggers for saving (save_now) and exiting (exit_now) !
!     operations, allowing user to force saving or gracefult exit from the     !
!     Process at any time during the simulation.                               !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Results_Type) :: Results    !! parent class
  integer             :: n_dom      !! number of domains
  type(Turb_Type)     :: Turb(MD)   !! turbulence modelling variables
  type(Field_Type)    :: Flow(MD)   !! flow fields
  type(Grid_Type)     :: Grid(MD)   !! computational grids
!----------------------------------[Locals]------------------------------------!
  integer :: d                   ! domain counter
  logical :: exit_now, save_now
!==============================================================================!

  ! Get out of here if you don't want to save initial conditions all right
  if( .not. Results % initial .and.  &
      Time % First_Dt() .eq. Time % Curr_Dt() ) then
    return
  end if

  inquire(file='exit_now', exist=exit_now)
  inquire(file='save_now', exist=save_now)

  ! Is it time to save results for post-processing?
  if(Time % Curr_Dt() .eq. Time % Last_Dt() .or.  &
     save_now                               .or.  &
     exit_now                               .or.  &
     Results % Time_To_Save_Results()       .or.  &
     Info % Time_To_Exit()) then

    do d = 1, n_dom
      call Control % Switch_To_Domain(d)
      call Results % Save_Vtu_Fields(Turb(d), Flow(d), Grid(d),      &
                                     plot_inside=.true., domain=d)
      call Results % Save_Vtu_Fields(Turb(d), Flow(d), Grid(d),      &
                                     plot_inside=.false., domain=d)

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
