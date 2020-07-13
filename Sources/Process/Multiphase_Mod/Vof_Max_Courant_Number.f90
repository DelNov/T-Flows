!==============================================================================!
  subroutine Vof_Max_Courant_Number(mult, dt, c_d, interf, courant_max)
!------------------------------------------------------------------------------!
!   Computes the Maximum Courant Number at cells. The argument interf helps    !
!   selecting if calculation will be performed close the interface, which      !
!   in turn will modify the time step if necessary and perform a simple        !
!   stepping approach. If interf = 0, then calculation is made everywhere      !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Multiphase_Type), target :: mult
  real                          :: dt     ! time step
  real                          :: c_d(-mult % pnt_grid % n_bnd_cells:  &
                                        mult % pnt_grid % n_cells) ! Courant n.
  integer                       :: interf
  real                          :: courant_max
!--------------------------------[Locals]--------------------------------------!
  type(Field_Type),     pointer :: flow
  type(Grid_Type),      pointer :: grid
  type(Var_Type),       pointer :: vof
  type(Face_Type),      pointer :: m_flux
  integer                       :: c, c1, c2, s
  real                          :: vof_dist
!==============================================================================!

  ! Take aliases
  flow     => mult % pnt_flow
  m_flux   => flow % m_flux
  grid     => flow % pnt_grid
  vof      => mult % vof

  courant_max = -HUGE

  ! Initialize
  c_d(-mult % pnt_grid % n_bnd_cells:mult % pnt_grid % n_cells) = 0.0

  if (interf == 1) then

    ! At boundaries
    do s = 1, grid % n_bnd_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)

      vof_dist = min(max(vof % n(c1), 0.0),1.0)

      vof_dist = (1.0 - vof_dist) * (1.0 - vof_dist)            &
                                  * vof_dist * vof_dist * 16.0

      c_d(c1) = c_d(c1) + vof_dist * max(- m_flux % n(s)              &
                                         / flow % density_f(s)        &
                                         * dt / grid % vol(c1), 0.0)
    end do

    ! Interior cells
    do s = grid % n_bnd_faces + 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)

      vof_dist = min(max(vof % n(c1), 0.0),1.0)

      vof_dist = (1.0 - vof_dist) * (1.0 - vof_dist)            &
                                  * vof_dist * vof_dist * 16.0

      c_d(c1) = c_d(c1) + vof_dist * max(- m_flux % n(s)              &
                                         / flow % density_f(s)        &
                                         * dt / grid % vol(c1), 0.0)

      vof_dist = min(max(vof % n(c2), 0.0),1.0)

      vof_dist = (1.0 - vof_dist) * (1.0 - vof_dist)            &
                                  * vof_dist * vof_dist * 16.0

      c_d(c2) = c_d(c2) + vof_dist * max( m_flux % n(s)              &
                                        / flow % density_f(s)        &
                                        * dt / grid % vol(c2), 0.0)

    end do

    !if (mult % phase_Change) then
    !  do c = 1, grid % n_cells
    !    vof_dist = min(max(vof % n(c1), 0.0),1.0)

    !    vof_dist = (1.0 - vof_dist) * (1.0 - vof_dist)            &
    !                                * vof_dist * vof_dist * 16.0

    !    c_d(c) = c_d(c) + vof_dist * mult % flux_rate(c)    &
    !                               / flow % density_f(s) * dt
    !  end do
    !end if

    call Grid_Mod_Exchange_Cells_Real(grid, c_d)

    do c = 1, grid % n_cells
      courant_max = max(c_d(c), courant_max)
    end do
    call Comm_Mod_Global_Max_Real(courant_max)

  else

    ! At boundaries
    do s = 1, grid % n_bnd_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)

      c_d(c1) = c_d(c1) + max(- m_flux % n(s) / flow % density_f(s)  &
                                        * dt / grid % vol(c1), 0.0)
    end do

    ! At interior faces
    do s = grid % n_bnd_faces + 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)

      c_d(c1) = c_d(c1) + max(- m_flux % n(s) / flow % density_f(s)  &
                                        * dt / grid % vol(c1), 0.0)

      c_d(c2) = c_d(c2) + max( m_flux % n(s) / flow % density_f(s)   &
                                        * dt / grid % vol(c2), 0.0)

    end do

    !if (mult % phase_Change) then
    !  do c = 1, grid % n_cells
    !    c_d(c) = c_d(c) + mult % flux_rate(c) / flow % density_f(s) * dt
    !  end do
    !end if

    call Grid_Mod_Exchange_Cells_Real(grid, c_d)
  end if

  end subroutine
