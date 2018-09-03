!==============================================================================!
  subroutine Calculate_Vis_T_Spalart_Allmaras(grid) 
!------------------------------------------------------------------------------!
!   Computes the turbulent viscosity for RANS models.                          !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Flow_Mod
  use Comm_Mod
  use Les_Mod
  use Rans_Mod
  use Grid_Mod
  use Control_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer         :: c
  real            :: x_rat, f_v1
!==============================================================================!

  if(turbulence_model .eq. DES_SPALART) then
    do c = 1, grid % n_cells
      x_rat    = vis % n(c)/viscosity
      f_v1     = x_rat**3/(x_rat**3 + c_v1**3)
      vis_t(c) = density * f_v1 * vis % n(c)
    end do
  end if

  if(turbulence_model .eq. SPALART_ALLMARAS) then
    do c = 1, grid % n_cells
      x_rat    = vis % n(c)/viscosity
      f_v1     = x_rat**3/(x_rat**3 + c_v1**3)
      vis_t(c) = density * f_v1 * vis % n(c)
    end do
  end if

  call Comm_Mod_Exchange_Real(grid, vis_t)  

  end subroutine
