!==============================================================================!
  subroutine Multiphase_Mod_Update_Physical_Properties(mult)
!------------------------------------------------------------------------------!
!   Update physical properties based on vof variable                           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Multiphase_Type), target :: mult
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: grid
  type(Field_Type), pointer :: flow
  type(Var_Type),   pointer :: vof
  real,             pointer :: vof_f(:)
  integer                   :: c, s
!==============================================================================!

  ! Take aliases
  flow  => mult % pnt_flow
  grid  => flow % pnt_grid
  vof   => mult % vof
  vof_f => mult % vof_f

  ! Density and viscosity in cells
  do c = 1, grid % n_cells
    density(c)   = vof % n(c)         * phase_dens(1)      &
                 + (1.0 - vof % n(c)) * phase_dens(2)
    viscosity(c) = vof % n(c)         * phase_visc(1)      &
                 + (1.0 - vof % n(c)) * phase_visc(2)
  end do
  call Comm_Mod_Exchange_Real(grid, density)
  call Comm_Mod_Exchange_Real(grid, viscosity)

  ! Density at faces
  do s = 1, grid % n_faces
    dens_face(s) =        vof_f(s)  * phase_dens(1)   &
                 + (1.0 - vof_f(s)) * phase_dens(2)
  end do

  end subroutine
