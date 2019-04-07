!==============================================================================!
  subroutine Calculate_Sgs_Hybrid(grid)
!------------------------------------------------------------------------------!
!   Calculates SGS stresses and turbulent viscosity for 'LES'.                 !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
  use Field_Mod
  use Comm_Mod
  use Turbulence_Mod
  use Grid_Mod
!------------------------------------------------------------------------------!
!   near(c) is the number of corresponding cell on the nearest wall.
!   In case that, in parallel executions, the subdomain does not have 
!   any nearwall cells, the near(c) is zero.
!   near(c) is calculated in NearWallCells.f90, only ones in the beginig
!   of a simulation.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer :: c
  real    :: lf
!==============================================================================!

  do c = 1, grid % n_cells
    lf = grid % vol(c) ** ONE_THIRD
    vis_t_sgs(c) = density                &
                 * (lf*lf)                &          ! delta^2 
                 * c_dyn(c)               &          ! c_dynamic   
                 * shear(c)
  end do

  call Comm_Mod_Exchange_Real(grid, vis_t_sgs)

  end subroutine
