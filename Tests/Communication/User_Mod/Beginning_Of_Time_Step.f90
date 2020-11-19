!==============================================================================!
  subroutine User_Mod_Beginning_Of_Time_Step(flow, turb, mult, swarm, n, time)
!------------------------------------------------------------------------------!
!   This function is called at the beginning of time step.                     !
!------------------------------------------------------------------------------!
  use Work_Mod, only: var => r_cell_01
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),      target :: flow
  type(Turb_Type),       target :: turb
  type(Multiphase_Type), target :: mult
  type(Swarm_Type),      target :: swarm
  integer, intent(in)           :: n     ! time step
  real,    intent(in)           :: time  ! physical time
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: grid
  type(Var_Type),  pointer :: u, v, w, t, phi, vof
  integer                  :: c, nb, nc
!==============================================================================!

  ! Take aliases
  grid => flow % pnt_grid
  nb   =  grid % n_bnd_cells
  nc   =  grid % n_cells

  do c = 1, grid % n_cells - grid % comm % n_buff_cells
    var(c) = grid % comm % cell_proc(c)
  end do

  call Grid_Mod_Save_Debug_Vtu(grid, 'var_before',                     &
                               scalar_cell = var(-nb:nc), scalar_name = 'var')

  call Grid_Mod_Exchange_Cells_Real(grid, var)

  call Grid_Mod_Save_Debug_Vtu(grid, 'var_after',                     &
                               scalar_cell = var(-nb:nc), scalar_name = 'var')

  call Comm_Mod_End
  stop

  end subroutine
