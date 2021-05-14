include '../User_Mod/Vof_Initialization_Ellipsoid.f90'

!==============================================================================!
  subroutine User_Mod_Initialize_Variables(flow, turb, Vof, swarm, sol)
!------------------------------------------------------------------------------!
!   User initialization of dependent variables.                                !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),  target :: flow
  type(Turb_Type),   target :: turb
  type(Vof_Type),    target :: Vof
  type(Swarm_Type),  target :: swarm
  type(Solver_Type), target :: sol
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: grid
  type(Var_Type),   pointer :: fun
  real,             pointer :: dt
  integer                   :: c
!==============================================================================!

  ! Take aliases
  grid => flow % pnt_grid
  fun  => Vof % fun
  dt   => flow % dt

  ! Initialize the whole domain as 0.0
  do c = 1, grid % n_cells
    fun % n(c) = 0.0
  end do

  ! Ellipsoid:
  call Vof_Initialization_Ellipsoid(Vof)

  call Grid_Mod_Exchange_Cells_Real(grid, fun % n)

  ! Old value
  fun % o(:) = fun % n(:)

  end subroutine
