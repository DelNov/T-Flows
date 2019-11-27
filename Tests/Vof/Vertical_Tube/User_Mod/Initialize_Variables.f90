include '../User_Mod/Vof_Initialization_Cylinder.f90'
include '../User_Mod/Vof_Initialization_Ellipsoid.f90'
include '../User_Mod/Vof_Initialization_Plane.f90'

!==============================================================================!
  subroutine User_Mod_Initialize_Variables(flow, turb, mult, swarm)
!------------------------------------------------------------------------------!
!   Case-dependent initialization of VOF variable.                             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),      target :: flow
  type(Turb_Type),       target :: turb
  type(Multiphase_Type), target :: mult
  type(Swarm_Type),      target :: swarm
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: grid
  type(Var_Type),   pointer :: vof
  real,             pointer :: dt
  real,             pointer :: vof_f(:)
  integer                   :: c, c1, c2, s
  real                      :: fs
!==============================================================================!

  ! Take aliases
  grid  => flow % pnt_grid
  vof   => mult % vof
  dt    => flow % dt
  vof_f => mult % vof_f

  ! Initialize the whole domain as 0.0
  do c = 1, grid % n_cells
    vof % n(c) = 0.0
  end do

  ! Under a Plane:
  ! call Vof_Initialization_Plane(mult)
  ! Ellipsoid:
  call Vof_Initialization_Ellipsoid(mult)
  ! Cylinder:
  ! call Vof_Initialization_Cylinder(mult)

  call Comm_Mod_Exchange_Real(grid, vof % n)

  ! Old value
  vof % o(:) = vof % n(:)

  end subroutine
