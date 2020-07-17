!==============================================================================!
  subroutine Multiphase_Mod_Update_Physical_Properties(mult, backup)
!------------------------------------------------------------------------------!
!   Update physical properties based on volume fraction variable               !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Multiphase_Type), target :: mult
  logical                       :: backup
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: grid
  type(Field_Type), pointer :: flow
  type(Var_Type),   pointer :: vof
  real, contiguous, pointer :: vof_f(:)
  integer                   :: c, s, c1, c2
  real                      :: fs
!==============================================================================!

  ! Take aliases
  flow  => mult % pnt_flow
  grid  => flow % pnt_grid
  vof   => mult % vof
  vof_f => mult % vof_f

  ! Density and viscosity in cells
  do c = 1, grid % n_cells
    flow % density(c)      = vof % n(c)         * mult % phase_dens(1)      &
                           + (1.0 - vof % n(c)) * mult % phase_dens(2)
    flow % viscosity(c)    = vof % n(c)         * mult % phase_visc(1)      &
                           + (1.0 - vof % n(c)) * mult % phase_visc(2)
    flow % capacity(c)     = vof % n(c)         * mult % phase_capa(1)      &
                           + (1.0 - vof % n(c)) * mult % phase_capa(2)
    flow % conductivity(c) = vof % n(c)         * mult % phase_cond(1)      &
                           + (1.0 - vof % n(c)) * mult % phase_cond(2)
  end do

  ! At boundaries

  do s = 1, grid % n_bnd_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    fs = grid % f(s)
    flow % viscosity(c2) = flow % viscosity(c1)
    flow % conductivity(c2) = flow % conductivity(c1)
  end do

  call Grid_Mod_Exchange_Cells_Real(grid, flow % density)
  call Grid_Mod_Exchange_Cells_Real(grid, flow % viscosity)

  if (backup .eqv. .false.) then
    ! Density at faces
    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)
      fs = grid % f(s)
      flow % density_f(s) =        vof_f(s)  * mult % phase_dens(1)   &
                          + (1.0 - vof_f(s)) * mult % phase_dens(2)
    end do

  end if

  end subroutine
