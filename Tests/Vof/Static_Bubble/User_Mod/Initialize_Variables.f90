include '../User_Mod/Vof_Initialization_Ellipsoid.f90'
include '../User_Mod/Vof_Interface_Ellipsoid.f90'

!==============================================================================!
  subroutine User_Mod_Initialize_Variables(Flow, turb, Vof, swarm, sol)
!------------------------------------------------------------------------------!
!   Case-dependent initialization of VOF variable.                             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),  target :: Flow
  type(Turb_Type),   target :: turb
  type(Vof_Type),    target :: Vof
  type(Swarm_Type),  target :: swarm
  type(Solver_Type), target :: sol
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: grid
  type(Var_Type),   pointer :: fun
  real,             pointer :: dt
  integer                   :: c, c1, c2, s
!==============================================================================!

  ! Take aliases
  grid => Flow % pnt_grid
  fun  => Vof % fun
  dt   => Flow % dt

  !---------------------------------!
  !   Initialize the VOF function   !
  !---------------------------------!

  ! Initialize the whole domain as 0.0
  fun % n(:) = 0.0

  ! Ellipsoid:
  call Vof_Initialization_Ellipsoid(Vof)

  ! Naive way to update bounary values
  do s = 1, grid % n_bnd_cells
    c1 = grid % faces_c(1, s)
    c2 = grid % faces_c(2, s)
    if(c2 < 0) then
      fun % n(c2) = fun % n(c1)
    end if
  end do

  ! Update buffer values
  call Grid_Mod_Exchange_Cells_Real(grid, fun % n)

  ! Set old values to be the same as new ones
  fun % o(:) = fun % n(:)

  !--------------------------------!
  !   Initialize front if needed   !
  !--------------------------------!
  ! if(Vof % track_front) then
  !   call Surf_Mod_Allocate(Vof % surf, Flow)
  !   call Surf_Mod_Place_At_Var_Value(Vof % surf,  &
  !                                    Vof % fun,   &
  !                                    sol,          &
  !                                    0.5,          &
  !                                    .false.)  ! don't print messages
  !   call Surf_Mod_Calculate_Curvatures_From_Elems(Vof % surf)
  !   call Surf_Mod_Compute_Distance_Function_And_Vof(Vof % surf,       &
  !                                                   Vof % dist_func,  &
  !                                                   Vof % fun)
  ! end if

  end subroutine
