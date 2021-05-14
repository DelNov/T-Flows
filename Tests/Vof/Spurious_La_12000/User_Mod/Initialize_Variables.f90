include '../User_Mod/Vof_Initialization_Cylinder.f90'
include '../User_Mod/Vof_Interface_Cylinder.f90'
include '../User_Mod/Vof_Area_Square_Circle.f90'
include '../User_Mod/Vof_Exact_Cylinder.f90'

!==============================================================================!
  subroutine User_Mod_Initialize_Variables(flow, turb, Vof, swarm, sol)
!------------------------------------------------------------------------------!
!   Case-dependent initialization of VOF variable.                             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),  target :: flow
  type(Turb_Type),   target :: turb
  type(Vof_Type),    target :: Vof
  type(Swarm_Type),  target :: swarm
  type(Solver_Type), target :: sol
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: grid
  type(Var_Type),  pointer :: fun
  real,            pointer :: dt
  integer                  :: c, c1, c2, s
!==============================================================================!

  ! Take aliases
  grid => flow % pnt_grid
  fun  => Vof % fun
  dt   => flow % dt

  ! Initialize the whole domain as 0.0
  do c = 1, grid % n_cells
    fun % n(c) = 0.0
  end do

  ! Cylinder:
  call Vof_Initialization_Cylinder(Vof)

  call Grid_Mod_Exchange_Cells_Real(grid, fun % n)

  ! Old value
  fun % o(:) = fun % n(:)

  ! At boundaries
  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    if(c2 < 0) then
      if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. OUTFLOW) then
        fun % n(c2) = fun % n(c1)
      else if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. INFLOW) then
      else
        fun % n(c2) = fun % n(c1)
      end if
    end if
  end do

  end subroutine
