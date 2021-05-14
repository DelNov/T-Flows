include '../User_Mod/Check_Inside_Box.f90'
include '../User_Mod/Vof_Initialization_Box.f90'
include '../User_Mod/Vof_Interface_Box.f90'

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
  type(Grid_Type),  pointer :: grid
  type(Var_Type),   pointer :: fun, t
  real,             pointer :: dt
  integer                   :: c, c1, c2, s
!==============================================================================!

  ! Take aliases
  grid => flow % pnt_grid
  t    => flow % t
  fun  => Vof % fun
  dt   => flow % dt

  !---------------------------------!
  !   Initialize the VOF function   !
  !---------------------------------!

  ! Initialize the whole domain as 0.0
  fun % n(:) = 0.0

  ! Box
  call Vof_Initialization_Box(Vof)

  ! Naive way to update bounary values
  do s = 1, grid % n_bnd_cells
    c1 = grid % faces_c(1, s)
    c2 = grid % faces_c(2, s)
    if(c2 < 0) then
      fun % n(c2) = fun % n(c1)
    end if
  end do

  ! Initialize front
  if(Vof % track_front) then
    call Front_Mod_Place_At_Var_Value(Vof % front,  &
                                      Vof % fun,    &
                                      sol,           &
                                      0.5,           &
                                      .true.)  ! don't print messages
    call Front_Mod_Calculate_Curvatures_From_Elems(Vof % front)
    call Front_Mod_Print_Statistics               (Vof % front)
    call Front_Mod_Save(Vof % front, 0)
  end if

  ! Initialize velocities (depends on phase definition)
  do c = 1, grid % n_cells

    ! Inside the water
    if(fun % n(c) .gt. 0.99) then
!     flow % u % n(c) =   0.05
!     flow % u % o(c) =   0.05
      if(flow % heat_transfer) then
        flow % t % n(c) = 100.0
        flow % t % o(c) = 100.0
      end if

    ! Inside the vapor
    else
      flow % u % n(c) =   0.0
      flow % u % o(c) =   0.0

      if(flow % heat_transfer) then
        ! small: flow % t % n(c) = 110.0 - grid % xc(c)/5.0e-5 * 10.0
        ! small: flow % t % o(c) = 110.0 - grid % xc(c)/5.0e-5 * 100.0
        ! MINI_1 flow % t % n(c) = 110.0 - grid % xc(c) * 100.0
        ! MINI_1 flow % t % o(c) = 110.0 - grid % xc(c) * 100.0
        ! MINI_2 flow % t % n(c) = 110.0 - grid % xc(c) * 400.0
        ! MINI_2 flow % t % o(c) = 110.0 - grid % xc(c) * 400.0
        ! MINI_3 flow % t % n(c) = 110.0 - grid % xc(c) * 4000.0
        ! MINI_3 flow % t % o(c) = 110.0 - grid % xc(c) * 4000.0
        ! MINI_4 flow % t % n(c) = 110.0 - grid % xc(c) * 10000.0
        ! MINI_4 flow % t % o(c) = 110.0 - grid % xc(c) * 10000.0
        flow % t % n(c) = 110.0 - grid % xc(c) * 20000.0
        flow % t % o(c) = 110.0 - grid % xc(c) * 20000.0
      end if
    end if
  end do

  ! Update buffer values
  call Grid_Mod_Exchange_Cells_Real(grid, fun % n)

  ! Set old values to be the same as new ones
  fun % o(:) = fun % n(:)

  end subroutine
