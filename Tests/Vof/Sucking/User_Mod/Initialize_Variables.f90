include '../User_Mod/Check_Inside_Box.f90'
include '../User_Mod/Vof_Initialization_Box.f90'
include '../User_Mod/Vof_Interface_Box.f90'

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
  type(Grid_Type),  pointer :: Grid
  type(Var_Type),   pointer :: fun, t
  real,             pointer :: dt
  real                      :: xc, tc
  integer                   :: c, c1, c2, s, l, fu
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
    call Vof % Smooth_For_Curvature_Csf()
    call Vof % Front % Place_Front_At_Value(Vof % fun,     &
                                            Vof % smooth,  &
                                            .true.)  ! don't print messages
    call Vof % Front % Print_Front_Statistics()
  end if

  ! Initialize temperatures from file

  call File % Open_For_Reading_Ascii("temperature_distribution_000.dat", fu)
  do l = 1, 100
    read(fu, *) xc, tc
    do c = 1, Grid % n_cells
      if(Math % Approx_Real(Grid % xc(c), xc)) then
        Flow % t % n(c) = tc
        Flow % t % o(c) = tc
      end if
    end do
  end do
  close(fu)

  ! Update buffer values
  call Grid % Exchange_Cells_Real(fun % n)

  ! Set old values to be the same as new ones
  fun % o (:) = fun % n(:)
  fun % oo(:) = fun % o(:)

  end subroutine
