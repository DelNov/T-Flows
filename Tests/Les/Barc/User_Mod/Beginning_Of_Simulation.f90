!==============================================================================!
  subroutine User_Mod_Beginning_Of_Simulation(Flow, Turb, Vof, Swarm,  &
                                              curr_dt, time)
!------------------------------------------------------------------------------!
!   This function is called at the beginning of simulation.                    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),    target :: Flow
  type(Turb_Type),     target :: Turb
  type(Vof_Type),      target :: Vof
  type(Swarm_Type),    target :: Swarm
  integer, intent(in)         :: curr_dt  ! time step
  real,    intent(in)         :: time     ! physical time
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
!==============================================================================!

  ! Make sure you compiled it
  if(this_proc < 2) print *, '# In User_Mod_Beginning_Of_Simulation'

  ! Take alias
  Grid => Flow % pnt_grid

  !--------------------------------------------------!
  !   Check the computed cells' moments of inertia   !
  !--------------------------------------------------!
  call Grid % Save_Debug_Vtu('cell-inertia',               &
                             tensor_cell = (/Grid % xx,    &
                                             Grid % yy,    &
                                             Grid % zz,    &
                                             Grid % xy,    &
                                             Grid % xz,    &
                                             Grid % yz/),  &
                             tensor_name = 'Cell Inertia')

  call Comm_Mod_End
  stop

  end subroutine
