!==============================================================================!
  subroutine User_Mod_Insert_Particles(Flow, Turb, Vof, Swarm, curr_dt, time)
!------------------------------------------------------------------------------!
!   Insert particles for Lagrangian particle tracking                          !
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
  type(Var_Type),  pointer :: u, v, w, t, phi
!==============================================================================!

  ! Take aliases
  Grid => Flow % pnt_grid

  ! Remove the following six lines in real Lagrangian tracking simulations
  if(this_proc < 2) then
    print *, '#======================================================'     //  &
             '======================================================='
    print *, '# WARNING: You are running a Lagrangian tracking simulation' //  &
             'with the default version of Insert_Particles.'
    print *, '# You have probably forgotten to compile Process with '      //  &
             'DIR_CASE=<full_or_relative_path_to_case> directive.'
    print *, '#------------------------------------------------------'     //  &
             '-------------------------------------------------------'
  end if

  end subroutine
