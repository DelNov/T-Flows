include '../User_Mod/Vof_Initialization_Cylinder.f90'
include '../User_Mod/Vof_Initialization_Ellipsoid.f90'
include '../User_Mod/Vof_Initialization_Plane.f90'
include '../User_Mod/Vof_Init_Random_Seed.f90'
include '../User_Mod/Vof_Interface_Cylinder.f90'
include '../User_Mod/Vof_Interface_Ellipsoid.f90'
include '../User_Mod/Vof_Interface_Plane.f90'
include '../User_Mod/Check_Inside_Cell.f90'
include '../User_Mod/Convert_Problem_Name_To_Integer.f90'
include '../User_Mod/Vof_Area_Square_Circle.f90'
include '../User_Mod/Vof_Exact_Cylinder.f90'
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
  !call Vof_Initialization_Cylinder(mult)

  call Grid_Mod_Exchange_Cells_Real(grid, vof % n)
  ! call Grid_Mod_Exchange_Real(grid, vof % n)

  ! Old value
  vof % o(:) = vof % n(:)

  do s = 1, grid % n_bnd_cells
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    fs = grid % f(s)
    if(Grid_Mod_Bnd_Cond_Type(grid,c2) .ne. INFLOW) then
      vof_f(s) = vof % n(c1)
    else
      vof_f(s) = vof % n(c2)
    end if
  end do

  do s = grid % n_bnd_cells + 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    fs = grid % f(s)
    vof_f(s) = fs * vof % n(c1) + (1.0 - fs) * vof % n(c2)
  end do
  end subroutine
