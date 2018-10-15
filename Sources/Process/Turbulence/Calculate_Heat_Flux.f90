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
  use Flow_Mod 
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
  real    :: beta, pr
!==============================================================================!

  call Grad_Mod_For_Phi(grid, t % n, 1, t_x, .true.)
  call Grad_Mod_For_Phi(grid, t % n, 2, t_y, .true.)
  call Grad_Mod_For_Phi(grid, t % n, 3, t_z, .true.)

  !-----------------------------------------!
  !   Compute the sources in the interior   !
  !-----------------------------------------!
  call Control_Mod_Turbulent_Prandtl_Number(pr_t)
  pr   = 0.71  ! bad, hard coded
  beta = 1.0

  if(turbulence_model .eq. K_EPS        .or.  &
     turbulence_model .eq. K_EPS_ZETA_F .or.  &
     turbulence_model .eq. DES_SPALART) then

    do c = 1, grid % n_cells
      pr_t = max(Turbulent_Prandtl_Number(grid, c), TINY)
      ut % n(c) = -vis_t(c) / pr_t * t_x(c)
      vt % n(c) = -vis_t(c) / pr_t * t_y(c)
      wt % n(c) = -vis_t(c) / pr_t * t_z(c)
    end do

  else if(turbulence_model .eq. RSM_MANCEAU_HANJALIC .or.  &
          turbulence_model .eq. RSM_HANJALIC_JAKIRLIC) then
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
  end if

  if(buoyancy) then
    do c = 1, grid % n_cells
      ut % n(c) = min( 0.01 * t_ref, ut % n(c))
      ut % n(c) = max(-0.01 * t_ref, ut % n(c))
      vt % n(c) = min( 0.01 * t_ref, vt % n(c))
      vt % n(c) = max(-0.01 * t_ref, vt % n(c))
      wt % n(c) = min( 0.01 * t_ref, wt % n(c))
      wt % n(c) = max(-0.01 * t_ref, wt % n(c))
      p_buoy(c) = -beta*(  grav_x * ut % n(c)  &
                         + grav_y * vt % n(c)  &
                         + grav_z * wt % n(c))
      p_buoy(c) = max(p_buoy(c),0.0)
    end do
  end if

  end subroutine
