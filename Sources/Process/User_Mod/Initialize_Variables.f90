!==============================================================================!
  subroutine User_Mod_Initialize_Variables(Flow, turb, Vof, swarm, Sol)
!------------------------------------------------------------------------------!
!   User initialization of dependent variables.                                !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),  target :: Flow
  type(Turb_Type),   target :: turb
  type(Vof_Type),    target :: Vof
  type(Swarm_Type),  target :: swarm
  type(Solver_Type), target :: Sol
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
  type(Var_Type),  pointer :: u, v, w, t, phi
!==============================================================================!

  ! Take aliases
  Grid => Flow % pnt_grid

  ! Remove the following 12 lines in real Lagrangian tracking simulations
  if(Vof % model .eq. VOLUME_OF_FLUID) then
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
