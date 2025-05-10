# ifdef __INTEL_COMPILER
#   include "User_Mod/Check_Inside_Box.f90"
#   include "User_Mod/Vof_Initialization_Box.f90"
#   include "User_Mod/Vof_Interface_Box.f90"
# else
#   include "Vof_Initialization_Box.f90"
#   include "Vof_Interface_Box.f90"
#   include "Check_Inside_Box.f90"
# endif

!==============================================================================!
  subroutine User_Mod_Initialize_Variables(Flow, Turb, Vof, Swarm, Sol)
!------------------------------------------------------------------------------!
!   Case-dependent initialization of VOF variable.                             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),  target :: Flow
  type(Turb_Type),   target :: Turb
  type(Vof_Type),    target :: Vof
  type(Swarm_Type),  target :: Swarm
  type(Solver_Type), target :: Sol
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: Grid
  type(Var_Type),   pointer :: fun, t
  real,             pointer :: dt
  integer                   :: c, c1, c2, s
!==============================================================================!

  ! Take aliases
  Grid => Flow % pnt_grid
  t    => Flow % t
  fun  => Vof % fun
  dt   => Flow % dt

  !---------------------------------!
  !   Initialize the VOF function   !
  !---------------------------------!

  ! Initialize the whole domain as 0.0
  fun % n(:) = 0.0

  ! Box
  call Vof_Initialization_Box(Vof)

  ! Naive way to update bounary values
  do s = 1, Grid % n_bnd_cells
    c1 = Grid % faces_c(1, s)
    c2 = Grid % faces_c(2, s)
    if(c2 < 0) then
      fun % n(c2) = fun % n(c1)
    end if
  end do

  ! Initialize front
  if(Vof % track_front) then
    call Vof % Smooth_Vof_And_Compute_Surface_Normals()
    call Vof % Front % Place_Front_At_Value(Vof % fun,     &
                                            .true.)  ! don't print messages
    call Vof % Front % Print_Front_Statistics()
  end if

  ! Initialize velocities (depends on phase definition)
  do c = Cells_In_Domain_And_Buffers()

    ! Inside the water
    if(fun % n(c) .gt. 0.99) then
!     Flow % u % n(c) =   0.05
!     Flow % u % o(c) =   0.05
      if(Flow % heat_transfer) then
        Flow % t % n(c) = 0.1
        Flow % t % o(c) = 0.1
      end if

    ! Inside the vapor
    else
      Flow % u % n(c) =   0.0
      Flow % u % o(c) =   0.0

      if(Flow % heat_transfer) then
        ! small: Flow % t % n(c) = 110.0 - Grid % xc(c)/5.0e-5 * 10.0
        ! small: Flow % t % o(c) = 110.0 - Grid % xc(c)/5.0e-5 * 100.0
        ! MINI_1 Flow % t % n(c) = 110.0 - Grid % xc(c) * 100.0
        ! MINI_1 Flow % t % o(c) = 110.0 - Grid % xc(c) * 100.0
        ! MINI_2 Flow % t % n(c) = 110.0 - Grid % xc(c) * 400.0
        ! MINI_2 Flow % t % o(c) = 110.0 - Grid % xc(c) * 400.0
        ! MINI_3 Flow % t % n(c) = 110.0 - Grid % xc(c) * 4000.0
        ! MINI_3 Flow % t % o(c) = 110.0 - Grid % xc(c) * 4000.0
        ! MINI_4 Flow % t % n(c) = 110.0 - Grid % xc(c) * 10000.0
        ! MINI_4 Flow % t % o(c) = 110.0 - Grid % xc(c) * 10000.0
        Flow % t % n(c) = 10.1 - Grid % xc(c) * 20000.0
        Flow % t % o(c) = 10.1 - Grid % xc(c) * 20000.0
      end if
    end if
  end do

  ! Update buffer values
  call Grid % Exchange_Cells_Real(fun % n)

  ! Set old values to be the same as new ones
  fun % o (:) = fun % n(:)
  fun % oo(:) = fun % o(:)

  end subroutine
