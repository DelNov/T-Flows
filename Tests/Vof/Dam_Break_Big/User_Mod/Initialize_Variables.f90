include '../User_Mod/Vof_Initialization_Cylinder.f90'
include '../User_Mod/Vof_Initialization_Ellipsoid.f90'
include '../User_Mod/Vof_Initialization_Plane.f90'
include '../User_Mod/Vof_Initialization_Box.f90'
include '../User_Mod/Vof_Init_Random_Seed.f90'
include '../User_Mod/Vof_Interface_Cylinder.f90'
include '../User_Mod/Vof_Interface_Ellipsoid.f90'
include '../User_Mod/Vof_Interface_Plane.f90'
include '../User_Mod/Vof_Interface_Box.f90'
include '../User_Mod/Check_Inside_Cell.f90'
include '../User_Mod/Check_Inside_Box.f90'
include '../User_Mod/Convert_Problem_Name_To_Integer.f90'
include '../User_Mod/Vof_Area_Square_Circle.f90'
include '../User_Mod/Vof_Exact_Cylinder.f90'
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
  integer                   :: c, c1, c2, s
!==============================================================================!

  ! Take aliases
  grid  => flow % pnt_grid
  vof   => mult % vof
  dt    => flow % dt

  !---------------------------------!
  !   Initialize the VOF function   !
  !---------------------------------!

  ! Initialize the whole domain as 0.0
  vof % n(:) = 0.0

  ! Under a Plane:
  ! call Vof_Initialization_Plane(mult)
  ! Ellipsoid:
  ! call Vof_Initialization_Ellipsoid(mult)
  ! Cylinder:
  ! call Vof_Initialization_Cylinder(mult)
  ! Box:
  call Vof_Initialization_Box(mult)

  ! Naive way to update bounary values
  do s = 1, grid % n_bnd_cells
    c1 = grid % faces_c(1, s)
    c2 = grid % faces_c(2, s)
    if(c2 < 0) then
      vof % n(c2) = vof % n(c1)
    end if
  end do

  ! Update buffer values
  call Grid_Mod_Exchange_Cells_Real(grid, vof % n)

  ! Set old values to be the same as new ones
  vof % o(:) = vof % n(:)

  end subroutine
