!==============================================================================!
  subroutine Turb_Mod_Vis_T_Hybrid_Les_Rans(turb)
!------------------------------------------------------------------------------!
!   Calculates SGS stresses and turbulent viscosity for 'LES'.                 !
!------------------------------------------------------------------------------!
!   near(c) is the number of corresponding cell on the nearest wall.
!   In case that, in parallel executions, the subdomain does not have 
!   any nearwall cells, the near(c) is zero.
!   near(c) is calculated in NearWallCells.f90, only ones in the beginig
!   of a simulation.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Turb_Type), target :: turb
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type), pointer :: flow
  type(Grid_Type),  pointer :: grid
  integer                   :: c
  real                      :: lf
!==============================================================================!

  ! Take aliases
  flow => turb % pnt_flow
  grid => flow % pnt_grid

  do c = 1, grid % n_cells
    lf = grid % vol(c) ** ONE_THIRD
    turb % vis_t_sgs(c) = density(c)       &
                        * (lf*lf)          &          ! delta^2 
                        * turb % c_dyn(c)  &          ! c_dynamic   
                        * flow % shear(c)
  end do

  call Comm_Mod_Exchange_Real(grid, turb % vis_t_sgs)

  end subroutine
