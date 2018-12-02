!==============================================================================!
  subroutine Calculate_Heat_Flux(grid)
!------------------------------------------------------------------------------!
!   Computes turbulent heat fluxes                                             !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
  use Control_Mod
  use Grid_Mod
  use Grad_Mod
  use Field_Mod
  use Rans_Mod
  use Work_Mod, only: t_x => r_cell_01,  &
                      t_y => r_cell_02,  &
                      t_z => r_cell_03
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!---------------------------------[Calling]------------------------------------!
  real :: Turbulent_Prandtl_Number
!-----------------------------------[Locals]-----------------------------------!
  integer :: c 
  real    :: beta
!==============================================================================!

  call Grad_Mod_For_Phi(grid, t % n, 1, t_x, .true.)
  call Grad_Mod_For_Phi(grid, t % n, 2, t_y, .true.)
  call Grad_Mod_For_Phi(grid, t % n, 3, t_z, .true.)

  !-----------------------------------------!
  !   Compute the sources in the interior   !
  !-----------------------------------------!
  call Control_Mod_Turbulent_Prandtl_Number(pr_t)

  if(turbulent_heat_flux_model .eq. SGDH) then

    do c = 1, grid % n_cells
      pr_t = max(Turbulent_Prandtl_Number(grid, c), TINY)
      ut % n(c) = -vis_t(c) / pr_t * t_x(c)
      vt % n(c) = -vis_t(c) / pr_t * t_y(c)
      wt % n(c) = -vis_t(c) / pr_t * t_z(c)
    end do

  else if(turbulent_heat_flux_model .eq. GGDH) then
   
    do c = 1, grid % n_cells

      ut % n(c) =  -0.22*t_scale(c) * (uu % n(c) * t_x(c) +  &
                                       uv % n(c) * t_y(c) +  &
                                       uw % n(c) * t_z(c))
      vt % n(c) =  -0.22*t_scale(c) * (uv % n(c) * t_x(c) +  &
                                       vv % n(c) * t_y(c) +  &
                                       vw % n(c) * t_z(c))
      wt % n(c) =  -0.22*t_scale(c) * (uw % n(c) * t_x(c) +  &
                                       vw % n(c) * t_y(c) +  &
                                       ww % n(c) * t_z(c))

    end do

  else if(turbulent_heat_flux_model .eq. AFM) then
  end if

  end subroutine
