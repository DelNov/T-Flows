!==============================================================================!
  subroutine Turb_Mod_Face_Vis(turb, vis_eff, s)
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
  type(Grid_Type),  pointer :: Grid
  type(Field_Type), pointer :: Flow
  integer                   :: c1, c2
!==============================================================================!

  ! Take alias
  Flow => turb % pnt_flow
  Grid => turb % pnt_grid

  c1 = Grid % faces_c(1,s)
  c2 = Grid % faces_c(2,s)

  ! Viscosity computed as a harmonic mean
  vis_eff = 2.0 / (    1.0 / Flow % viscosity(c1)  &
                     + 1.0 / Flow % viscosity(c2) )

  if(turb % model .ne. NO_TURBULENCE_MODEL .and.  &
     turb % model .ne. DNS                 .and.  &
     turb % model .ne. HYBRID_LES_RANS) then
    vis_eff = vis_eff + Grid % fw(s)  * turb % vis_t(c1)  &
                  +(1.0-Grid % fw(s)) * turb % vis_t(c2)
  end if

  if(turb % model .eq. HYBRID_LES_RANS) then
    vis_eff =      Grid % fw(s)  * turb % vis_t_eff(c1)   &
            + (1.0-Grid % fw(s)) * turb % vis_t_eff(c2) + vis_eff
  end if

  if(c2 < 0) then
    if( turb % model .eq. LES_SMAGORINSKY    .or.  &
        turb % model .eq. LES_DYNAMIC        .or.  &
        turb % model .eq. HYBRID_LES_PRANDTL .or.  &
        turb % model .eq. LES_WALE) then
      if(Grid % Bnd_Cond_Type(c2) .eq. WALL .or.  &
         Grid % Bnd_Cond_Type(c2) .eq. WALLFL) then
        vis_eff = turb % vis_w(c1)
      end if
    end if
  end if

  if( turb % model .eq. K_EPS_ZETA_F     .or.  &
      turb % model .eq. HYBRID_LES_RANS  .or.  &
      turb % model .eq. K_EPS) then
    if(c2 < 0) then
      if(Grid % Bnd_Cond_Type(c2) .eq. WALL .or.  &
         Grid % Bnd_Cond_Type(c2) .eq. WALLFL) then
        vis_eff = turb % vis_w(c1)
      end if
    end if
  end if

  end subroutine
