!==============================================================================!
  subroutine User_Mod_Initialize_Variables(Flow, Turb, Vof, Swarm, Sol)
!------------------------------------------------------------------------------!
!   User initialization of dependent variables.                                !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),  target :: Flow
  type(Turb_Type),   target :: Turb
  type(Vof_Type),    target :: Vof
  type(Swarm_Type),  target :: Swarm
  type(Solver_Type), target :: Sol
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
  type(Var_Type),  pointer :: u, v, w, t, phi
!==============================================================================!

  ! Take aliases
  Grid => Flow % pnt_grid

  ! Remove the following 12 lines in real Lagrangian tracking simulations
  if(Flow % with_interface) then
    if(this_proc < 2) then
      print *, '#======================================================' //  &
               '======================================================='
      print *, '# WARNING: You are running a Volume of Fluid simulation' //  &
               'with the default version of Initialize_Variables.'
      print *, '# You have probably forgotten to compile Process with '  //  &
               'DIR_CASE=<full_or_relative_path_to_case> directive.'
      print *, '#------------------------------------------------------' //  &
               '-------------------------------------------------------'
    end if
  end if

  end subroutine
