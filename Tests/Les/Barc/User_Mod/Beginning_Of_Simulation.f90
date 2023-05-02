!==============================================================================!
  subroutine User_Mod_Beginning_Of_Simulation(Flow, Turb, Vof, Swarm)
!------------------------------------------------------------------------------!
!   This function is called at the beginning of simulation.                    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),    target :: Flow
  type(Turb_Type),     target :: Turb
  type(Vof_Type),      target :: Vof
  type(Swarm_Type),    target :: Swarm
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
!==============================================================================!

  ! Make sure you compiled it
  if(First_Proc()) print *, '# In User_Mod_Beginning_Of_Simulation'

  ! Take alias
  Grid => Flow % pnt_grid

  !--------------------------------------------------!
  !   Check the computed cells' moments of inertia   !
  !--------------------------------------------------!
  call Grid % Save_Debug_Vtu('cell-inertia',                &
                             tensor_cell = (/Grid % ixx,    &
                                             Grid % iyy,    &
                                             Grid % izz,    &
                                             Grid % ixy,    &
                                             Grid % ixz,    &
                                             Grid % iyz/),  &
                             tensor_comp = 6,               &
                             tensor_name = 'Cell Inertia')

  call Global % End_Parallel
  stop

  end subroutine
