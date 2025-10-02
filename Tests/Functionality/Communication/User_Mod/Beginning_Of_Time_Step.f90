!==============================================================================!
  subroutine User_Mod_Beginning_Of_Time_Step(Flow, Turb, Vof, Swarm)
!------------------------------------------------------------------------------!
!   This function is called at the beginning of time step.                     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type), target :: Flow
  type(Turb_Type),  target :: Turb
  type(Vof_Type),   target :: Vof
  type(Swarm_Type), target :: Swarm
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: G
  type(Var_Type),   pointer :: u, v, w, t, phi
  integer                   :: c, nb, nc
  real, contiguous, pointer :: var(:)
!==============================================================================!

  call Work % Connect_Real_Cell(var)

  ! Take aliases
  G  => Flow % pnt_grid
  nb =  G % n_bnd_cells
  nc =  G % n_cells

  ! In x direction
  var(:) = 0.0
  do c = 1, G % n_cells - G % Comm % n_buff_cells
    var(c) = G % xc(c)
  end do
  call G % Save_Debug_Vtu('x1',scalar_cell=var(-nb:nc),scalar_name='x')
  call G % Exchange_Cells_Real(var)
  call G % Save_Debug_Vtu('x2',scalar_cell=var(-nb:nc),scalar_name='x')

  ! In y direction
  var(:) = 0.0
  do c = 1, G % n_cells - G % Comm % n_buff_cells
    var(c) = G % yc(c)
  end do
  call G % Save_Debug_Vtu('y1',scalar_cell=var(-nb:nc),scalar_name='y')
  call G % Exchange_Cells_Real(var)
  call G % Save_Debug_Vtu('y2',scalar_cell=var(-nb:nc),scalar_name='y')

  ! In z direction
  var(:) = 0.0
  do c = 1, G % n_cells - G % Comm % n_buff_cells
    var(c) = G % zc(c)
  end do
  call G % Save_Debug_Vtu('z1',scalar_cell=var(-nb:nc),scalar_name='z')
  call G % Exchange_Cells_Real(var)
  call G % Save_Debug_Vtu('z2',scalar_cell=var(-nb:nc),scalar_name='z')

  call Work % Disconnect_Real_Cell(var)

  call Global % End_Parallel
  stop

  end subroutine
