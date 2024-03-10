!==============================================================================!
  subroutine User_Mod_Insert_Particles(Flow, Turb, Vof, Swarm)
!------------------------------------------------------------------------------!
!   Insert particles for Lagrangian particle tracking                          !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type), target :: Flow
  type(Turb_Type),  target :: Turb
  type(Vof_Type),   target :: Vof
  type(Swarm_Type), target :: Swarm
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
  type(Var_Type),  pointer :: u, v, w, t, phi
!==============================================================================!

  ! Take aliases
  Grid => Flow % pnt_grid

  ! Remove the following six lines in real Lagrangian tracking simulations
  call Message % Error(72,                                               &
             'You are running a Lagrangian tracking simulation  '    //  &
             'with the default version of Insert_Particles. You '    //  &
             'have probably forgotten to compile Process with '      //  &
             'DIR_CASE=<full_or_relative_path_to_case> directive.',      &
             file=__FILE__, line=__LINE__, one_proc=.true.)

  end subroutine
