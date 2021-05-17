!==============================================================================!
  subroutine User_Mod_Beginning_Of_Time_Step(Flow, turb, Vof, swarm, n, time)
!------------------------------------------------------------------------------!
!   This function is called at the beginning of time step.                     !
!------------------------------------------------------------------------------!
  use Work_Mod, only: var => r_cell_01
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type), target :: Flow
  type(Turb_Type),  target :: turb
  type(Vof_Type),   target :: Vof
  type(Swarm_Type), target :: swarm
  integer, intent(in)      :: n     ! time step
  real,    intent(in)      :: time  ! physical time
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: g
  type(Var_Type),  pointer :: u, v, w, t, phi
  integer                  :: c, nb, nc
!==============================================================================!

  ! Take aliases
  g  => Flow % pnt_grid
  nb =  g % n_bnd_cells
  nc =  g % n_cells

  ! In x direction
  var(:) = 0.0
  do c = 1, g % n_cells - g % comm % n_buff_cells
    var(c) = g % xc(c)
  end do
  call Grid_Mod_Save_Debug_Vtu(g,'x1',scalar_cell=var(-nb:nc),scalar_name='x')
  call Grid_Mod_Exchange_Cells_Real(g, var)
  call Grid_Mod_Save_Debug_Vtu(g,'x2',scalar_cell=var(-nb:nc),scalar_name='x')

  ! In y direction
  var(:) = 0.0
  do c = 1, g % n_cells - g % comm % n_buff_cells
    var(c) = g % yc(c)
  end do
  call Grid_Mod_Save_Debug_Vtu(g,'y1',scalar_cell=var(-nb:nc),scalar_name='y')
  call Grid_Mod_Exchange_Cells_Real(g, var)
  call Grid_Mod_Save_Debug_Vtu(g,'y2',scalar_cell=var(-nb:nc),scalar_name='y')

  ! In z direction
  var(:) = 0.0
  do c = 1, g % n_cells - g % comm % n_buff_cells
    var(c) = g % zc(c)
  end do
  call Grid_Mod_Save_Debug_Vtu(g,'z1',scalar_cell=var(-nb:nc),scalar_name='z')
  call Grid_Mod_Exchange_Cells_Real(g, var)
  call Grid_Mod_Save_Debug_Vtu(g,'z2',scalar_cell=var(-nb:nc),scalar_name='z')

  call Comm_Mod_End
  stop

  end subroutine
