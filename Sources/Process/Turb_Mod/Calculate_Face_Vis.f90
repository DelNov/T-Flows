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
  type(Field_Type), pointer :: flow
  integer                   :: c1, c2
!==============================================================================!

  ! Take alias
  flow => turb % pnt_flow
  grid => turb % pnt_grid

  c1 = grid % faces_c(1,s)
  c2 = grid % faces_c(2,s)

  vis_eff =        grid % fw(s)  * flow % viscosity(c1)   &
          + (1.0 - grid % fw(s)) * flow % viscosity(c2)

  vis_eff = 2.0 / (    1.0 / flow % viscosity(c1)  &
                     + 1.0 / flow % viscosity(c2) )

  if(turb % model .ne. NO_TURBULENCE_MODEL .and.  &
     turb % model .ne. DNS                 .and.  &
     turb % model .ne. HYBRID_LES_RANS) then
    vis_eff = vis_eff + grid % fw(s)  * turb % vis_t(c1)  &
                  +(1.0-grid % fw(s)) * turb % vis_t(c2)
  end if

  if(turb % model .eq. HYBRID_LES_RANS) then
    vis_eff =      grid % fw(s)  * turb % vis_t_eff(c1)   &
            + (1.0-grid % fw(s)) * turb % vis_t_eff(c2) + vis_eff
  end if

  if(c2 < 0) then
    if( turb % model .eq. LES_SMAGORINSKY    .or.  &
        turb % model .eq. LES_DYNAMIC        .or.  &
        turb % model .eq. HYBRID_LES_PRANDTL .or.  &
        turb % model .eq. LES_WALE) then
      if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALL .or.  &
         Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALLFL) then
        vis_eff = turb % vis_w(c1)
      end if
    end if
  end if

  if( turb % model .eq. K_EPS_ZETA_F     .or.  &
      turb % model .eq. HYBRID_LES_RANS  .or.  &
      turb % model .eq. K_EPS) then 
    if(c2 < 0) then
      if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALL .or.  &
         Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALLFL) then
        vis_eff = turb % vis_w(c1)
      end if
    end if
  end if

  end subroutine
