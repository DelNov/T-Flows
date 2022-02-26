!==============================================================================!
  subroutine User_Mod_Beginning_Of_Time_Step(Flow, Turb, Vof, Swarm, n, time)
!------------------------------------------------------------------------------!
!   This function is called at the beginning of time step.                     !
!------------------------------------------------------------------------------!
  use Work_Mod, only: var => r_cell_01
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type), target :: Flow
  type(Turb_Type),  target :: Turb
  type(Vof_Type),   target :: Vof
  type(Swarm_Type), target :: Swarm
  integer, intent(in)      :: n     ! time step
  real,    intent(in)      :: time  ! physical time
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: G
  type(Var_Type),  pointer :: u, v, w, t, phi
  integer                  :: c, nb, nc
!==============================================================================!

  ! Take aliases
  G  => Flow % pnt_grid
  nb =  G % n_bnd_cells
  nc =  G % n_cells

  ! In x direction
  var(:) = 0.0
  do c = 1, G % n_cells - G % comm % n_buff_cells
    var(c) = G % xc(c)
  end do
  call G % Save_Debug_Vtu('x1',scalar_cell=var(-nb:nc),scalar_name='x')
  call G % Exchange_Cells_Real(var)
  call G % Save_Debug_Vtu('x2',scalar_cell=var(-nb:nc),scalar_name='x')

  ! In y direction
  var(:) = 0.0
  do c = 1, G % n_cells - G % comm % n_buff_cells
    var(c) = G % yc(c)
  end do
  call G % Save_Debug_Vtu('y1',scalar_cell=var(-nb:nc),scalar_name='y')
  call G % Exchange_Cells_Real(var)
  call G % Save_Debug_Vtu('y2',scalar_cell=var(-nb:nc),scalar_name='y')

  ! In z direction
  var(:) = 0.0
  do c = 1, G % n_cells - G % comm % n_buff_cells
    var(c) = G % zc(c)
  end do
  call G % Save_Debug_Vtu('z1',scalar_cell=var(-nb:nc),scalar_name='z')
  call G % Exchange_Cells_Real(var)
  call G % Save_Debug_Vtu('z2',scalar_cell=var(-nb:nc),scalar_name='z')

  call Comm_Mod_End
  stop

  end subroutine
