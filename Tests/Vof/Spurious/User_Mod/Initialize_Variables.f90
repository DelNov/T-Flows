include '../User_Mod/Vof_Initialization_Cylinder.f90'
include '../User_Mod/Vof_Exact_Cylinder.f90'
include '../User_Mod/Vof_Area_Square_Circle.f90'

!==============================================================================!
  subroutine User_Mod_Initialize_Variables(flow, turb, mult, swarm, sol)
!------------------------------------------------------------------------------!
!   Case-dependent initialization of VOF variable.                             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),      target :: flow
  type(Turb_Type),       target :: turb
  type(Multiphase_Type), target :: mult
  type(Swarm_Type),      target :: swarm
  type(Solver_Type),     target :: sol
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: grid
  type(Var_Type),   pointer :: vof
  real,             pointer :: dt
  real,             pointer :: vof_f(:)
  integer                   :: c, c1, c2, s, i_probe, c_inters, n, code
  integer                   :: i, j
  real                      :: fs, min_dist, dist, glo_dist, epsloc
!==============================================================================!

  ! Take aliases
  grid => flow % pnt_grid
  vof  => mult % vof
  dt   => flow % dt

  epsloc = epsilon(epsloc)

  ! Initialize the whole domain as 0.0
  do c = 1, grid % n_cells
    vof % n(c) = 0.0
  end do

  ! Under a Plane:
  ! call Vof_Initialization_Plane(mult)
  ! Ellipsoid:
  ! call Vof_Initialization_Ellipsoid(mult)
  ! Cylinder:
  call Vof_Initialization_Cylinder(mult)
  ! Box:
  ! call Vof_Initialization_Box(mult)

  call Grid_Mod_Exchange_Cells_Real(grid, vof % n)

  ! Old value
  vof % o(:) = vof % n(:)

  ! At faces

  ! At boundaries
  do s = 1, grid % n_bnd_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. OUTFLOW) then
      vof % n(c2) = vof % n(c1)
    else if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. INFLOW) then
    else
      vof % n(c2) = vof % n(c1)
    end if
  end do

  end subroutine
