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

  !----------------------------!
  !   Allocate user's fields   !
  !- - - - - - - - - - - - - - +-------------------!
  !   User-defined fields must be defined in the   !
  !   range: -Grid % n_bnd_cells:Grid % n_cells)   !
  !------------------------------------------------!
  !@ allocate(my_useful_field(-Grid % n_bnd_cells : Grid % n_cells))

  !---------------------------------------------------------------------------!
  !   Remove the following 12 lines in real Lagrangian tracking simulations   !
  !---------------------------------------------------------------------------!
  if(Flow % with_interface .and. .not. Vof % init_stl) then
    call Message % Error(96,                                                 &
                'You are running a Volume of Fluid simulation with the '  // &
                'default version of Initialize_Variables, and you did '   // &
                'not provide an STL file for initialization.  You might'  // &
                'have forgotten to compile Process with '                 // &
                'DIR_CASE=<full_or_relative_path_to_case> directive.',       &
                file=__FILE__, line=__LINE__, one_proc=.true.)
  end if

  end subroutine
