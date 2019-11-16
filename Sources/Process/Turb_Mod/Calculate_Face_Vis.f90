!==============================================================================!
  subroutine Turb_Mod_Calculate_Face_Vis(turb, vis_eff, s)
!------------------------------------------------------------------------------!
!   Computes turbulent viscosity on a cell face for all turbulence models.     !
!   It is called from Compute_Momentum, while discretizing diffusion terms.    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Turb_Type), target :: turb
  real                    :: vis_eff
  integer                 :: s
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: grid
  integer                   :: c1, c2
  real                      :: fs
!==============================================================================!

  ! Take alias
  grid => turb % pnt_grid

  c1 = grid % faces_c(1,s)
  c2 = grid % faces_c(2,s)
  fs = grid % f(s)

  if (c2 > 0) then
    vis_eff = fs * viscosity(c1) + (1.0 - fs) * viscosity(c2)
  else
    vis_eff = viscosity(c1)
  end if

  if(turbulence_model .ne. NONE .and.  &
     turbulence_model .ne. DNS) then
    vis_eff = vis_eff + grid % fw(s)  * turb % vis_t(c1)  &
                  +(1.0-grid % fw(s)) * turb % vis_t(c2)
  end if

  if(turbulence_model .eq. HYBRID_LES_RANS) then
    vis_eff =      grid % fw(s)  * turb % vis_t_eff(c1)   &
            + (1.0-grid % fw(s)) * turb % vis_t_eff(c2) + vis_eff
  end if

  if(c2 < 0) then
    if( turbulence_model .eq. LES_SMAGORINSKY    .or.  &
        turbulence_model .eq. LES_DYNAMIC        .or.  &
        turbulence_model .eq. HYBRID_LES_PRANDTL .or.  &
        turbulence_model .eq. LES_WALE) then
      if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALL .or.  &
         Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALLFL) then
        vis_eff = turb % vis_w(c1)
      end if
    end if
  end if

  if( turbulence_model .eq. K_EPS_ZETA_F     .or.  &
      turbulence_model .eq. HYBRID_LES_RANS  .or.  &
      turbulence_model .eq. K_EPS) then 
    if(c2 < 0) then
      if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALL .or.  &
         Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALLFL) then
        vis_eff = turb % vis_w(c1)
      end if
    end if
  end if

  end subroutine
