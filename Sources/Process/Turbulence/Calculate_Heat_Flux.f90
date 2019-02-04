!==============================================================================!
  subroutine Calculate_Heat_Flux(flow)
!------------------------------------------------------------------------------!
!   Computes turbulent heat fluxes                                             !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
  use Control_Mod
  use Grid_Mod,  only: Grid_Type
  use Grad_Mod
  use Field_Mod, only: Field_Type
  use Rans_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type), target :: flow
!---------------------------------[Calling]------------------------------------!
  real :: Turbulent_Prandtl_Number
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: grid
  type(Var_Type),  pointer :: t
  integer                  :: c
!==============================================================================!

  ! Take aliases
  grid => flow % pnt_grid
  t    => flow % t

  ! Check if these are already computed somewhere, ...
  ! ... maybe this call is not needed
  call Grad_Mod_Array(grid, t % n, t % x, t % y, t % z, .true.)

  !-----------------------------------------!
  !   Compute the sources in the interior   !
  !-----------------------------------------!
  call Control_Mod_Turbulent_Prandtl_Number(pr_t)

  if(turbulent_heat_flux_model .eq. SGDH) then

    do c = 1, grid % n_cells
      pr_t = max(Turbulent_Prandtl_Number(grid, c), TINY)
      ut % n(c) = -vis_t(c) / pr_t * t % x(c)
      vt % n(c) = -vis_t(c) / pr_t * t % y(c)
      wt % n(c) = -vis_t(c) / pr_t * t % z(c)
    end do

  else if(turbulent_heat_flux_model .eq. GGDH) then
   
    do c = 1, grid % n_cells
      ut % n(c) =  -c_theta*t_scale(c) * (uu % n(c) * t_x(c)  +  &
                                          uv % n(c) * t_y(c)  +  &
                                          uw % n(c) * t_z(c))
      vt % n(c) =  -c_theta*t_scale(c) * (uv % n(c) * t_x(c)  +  &
                                          vv % n(c) * t_y(c)  +  &
                                          vw % n(c) * t_z(c))
      wt % n(c) =  -c_theta*t_scale(c) * (uw % n(c) * t_x(c)  +  &
                                          vw % n(c) * t_y(c)  +  &
                                          ww % n(c) * t_z(c))
    end do

  else if(turbulent_heat_flux_model .eq. AFM) then
  end if

  end subroutine
