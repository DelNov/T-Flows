!#include "Pv_Sat.f90"

!==============================================================================!
  subroutine User_Mod_Initialize_Variables(Flow, Turb, Vof, Swarm, Sol)
!------------------------------------------------------------------------------!
!   User initialization of dependent variables.                                !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),  target :: Flow
  type(Turb_Type),   target :: Turb
  type(Vof_Type),    target :: Vof
  type(Swarm_Type),  target :: Swarm
  type(Solver_Type), target :: Sol
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
  type(Var_Type),  pointer :: t, scalar
  integer                  :: c, s, c1, c2
  real                     :: t_tmp, p_v, m_h2o, m_air
  real                     :: t_cold, t_hot, z_ext, a_salt
!==============================================================================!

  ! Take aliases
  Grid   => Flow % pnt_grid
  t      => Flow % t
  scalar => Flow % scalar(1)

  !--------------------------------!
  !   Hot channel initalizations   !
  !--------------------------------!
  if(Grid % name .eq. 'UPPER_DOM') then

    ! Hopefully read from membrane_input file
    ! one day until then it is hard-coded
    a_salt = 0.09259259 ! salt[g/l] / rho[kg/m3]

    ! Allocate membrane interface temperaturefrom previous iteration and flux
    allocate (t_int_mem_prev(-Grid % n_bnd_cells:Grid % n_cells))
    allocate (mem_j(-Grid % n_bnd_cells:Grid % n_cells))

    do c = Cells_In_Domain_And_Buffers()
      ! initialize membrane interface temperature
      t_int_mem_prev(c) = Flow % t % n(c)
      ! Salt concentration
      scalar % n(c)  = a_salt
    end do

    ! Salt concentration inlet boundary
    do s = 1, Grid % n_faces
      c2 = Grid % faces_c(2,s) !-> each face 2 cells, right cell
      if (c2 .lt. 0) then !at boundary

        ! Check boundary by name 
        if (Var_Mod_Bnd_Cond_Type(scalar,c2) .eq. INFLOW) then
          scalar % n(c2) = a_salt
        end if

      end if
    end do

  end if

  !------------------------------!
  !   Air-gap: initializations   !
  !------------------------------!
  if(Grid % name .eq. 'LOWER_DOM') then

    ! Allocate membrane interface temperature from previous iteration
    allocate (t_int_prev(-Grid % n_bnd_cells:Grid % n_cells))

    ! Hopefully read from membrane_input file
    ! one day until then it is hard-coded
    m_h2o  = 18e+3 ! kg/mol
    m_air  = 28e+3 ! kg/mol
    t_cold = 15.0  ! deg C
    t_hot  = 74.68  ! deg C

    ! Modify for multiple processors:
    z_ext  = maxval(Grid % zc(:))  &
           - minval(Grid % zc(:))  ! average height of the air-gap in m

    do c = Cells_In_Domain_And_Buffers()
      t_int_prev(c) = t_cold

      ! Temperature (linear profile between hot and cold T)
      Flow % t % n(c) = t_hot + (t_hot - t_cold) * Grid % zc(c) / z_ext

      ! Vapor content according to temperature
      call Pv_Sat(Flow % t % n(c), p_v)
      scalar % n(c)  = p_v * 1e-5 * m_h2o / m_air
    end do

    ! Initialize scalar boundary value at membrane
    ! important for the partial vapor pressure calculation
    ! in Interface_Exchange
    do s = 1, Grid % n_faces
      c1 = Grid % faces_c(1,s) !-> each face 2 cells, left cell
      c2 = Grid % faces_c(2,s) !-> each face 2 cells, right cell

      ! At boundary, check it by name
      if (c2 .lt. 0) then
        if (Var_Mod_Bnd_Cond_Name(scalar,c2) .eq. 'TOP_WALL') then
          scalar % n(c2) = scalar % n(c1)
        end if
      end if

    end do
  end if

  end subroutine
