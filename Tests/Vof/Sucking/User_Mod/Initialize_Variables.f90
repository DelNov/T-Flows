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
  real                      :: xc, tc
  integer                   :: c, c1, c2, s, l, fu, nc
  character(SL)             :: hash
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

  ! Initialize temperatures and VOF from file
  call File % Open_For_Reading_Ascii("temperature_distribution_000.exa", fu)
  read(fu, *) hash, nc
  do l = 1, nc
    read(fu, *) xc, tc
    do c = Cells_In_Domain_And_Buffers()
      if(Math % Approx_Real(Grid % xc(c), xc)) then
        Flow % t % n(c) = tc
        Flow % t % o(c) = tc
        if(tc < 10.0+MICRO) then
          fun % n(c) = 0.0
        else
          fun % n(c) = 1.0
        end if
      end if
    end do
  end do
  close(fu)

  ! Naive way to update bounary values
  do s = 1, Grid % n_bnd_cells
    c1 = Grid % faces_c(1, s)
    c2 = Grid % faces_c(2, s)
    if(c2 < 0) then
      fun % n(c2) = fun % n(c1)
    end if
  end do

  ! Update buffer values
  call Grid % Exchange_Cells_Real(fun % n)

  ! Set old values to be the same as new ones
  fun % o (:) = fun % n(:)
  fun % oo(:) = fun % o(:)

  ! Initialize front
  if(Vof % track_front) then
    call Vof % Smooth_Vof_And_Compute_Surface_Normals()
    call Vof % Front % Place_Front_At_Value(Vof % fun,     &
                                            .true.)  ! don't print messages
    call Vof % Front % Print_Front_Statistics()
  end if

  end subroutine
