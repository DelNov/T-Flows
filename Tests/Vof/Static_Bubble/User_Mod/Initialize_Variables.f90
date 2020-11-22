include '../User_Mod/Vof_Initialization_Ellipsoid.f90'
include '../User_Mod/Vof_Interface_Ellipsoid.f90'
include '../User_Mod/Check_Inside_Cell.f90'

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
  grid => flow % pnt_grid
  vof  => mult % vof
  dt   => flow % dt

  !---------------------------------!
  !   Initialize the VOF function   !
  !---------------------------------!

  ! Initialize the whole domain as 0.0
  vof % n(:) = 0.0

  ! Ellipsoid:
  call Vof_Initialization_Ellipsoid(mult)

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

  !--------------------------------!
  !   Initialize front if needed   !
  !--------------------------------!
  ! if(mult % track_front) then
  !   call Surf_Mod_Allocate(mult % surf, flow)
  !   call Surf_Mod_Place_At_Var_Value(mult % surf,  &
  !                                    mult % vof,   &
  !                                    sol,          &
  !                                    0.5,          &
  !                                    .false.)  ! don't print messages
  !   call Surf_Mod_Calculate_Curvatures_From_Elems(mult % surf)
  !   call Surf_Mod_Compute_Distance_Function_And_Vof(mult % surf,       &
  !                                                   mult % dist_func,  &
  !                                                   mult % vof)
  ! end if

  end subroutine
