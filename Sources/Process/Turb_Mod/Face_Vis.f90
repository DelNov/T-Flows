!==============================================================================!
  subroutine Face_Vis(Turb, vis_eff, s)
!------------------------------------------------------------------------------!
!   Computes turbulent viscosity on a cell face for all turbulence models.     !
!   It is called from Compute_Momentum, while discretizing diffusion terms.    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Turb_Type), target :: Turb
  real,        intent(out) :: vis_eff
  integer,     intent(in)  :: s
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: Grid
  type(Field_Type), pointer :: Flow
  integer                   :: c1, c2
!==============================================================================!

  ! Take alias
  Flow => Turb % pnt_flow
  Grid => Turb % pnt_grid

  c1 = Grid % faces_c(1,s)
  c2 = Grid % faces_c(2,s)

  ! Viscosity computed as a harmonic mean
  vis_eff = 2.0 / (    1.0 / Flow % viscosity(c1)  &
                     + 1.0 / Flow % viscosity(c2) )

  if(Turb % model .ne. NO_TURBULENCE_MODEL .and.  &
     Turb % model .ne. DNS                 .and.  &
     Turb % model .ne. HYBRID_LES_RANS) then
    vis_eff = vis_eff + Grid % fw(s)  * Turb % vis_t(c1)  &
                  +(1.0-Grid % fw(s)) * Turb % vis_t(c2)
  end if

  if(Turb % model .eq. HYBRID_LES_RANS) then
    vis_eff =      Grid % fw(s)  * Turb % vis_t_eff(c1)   &
            + (1.0-Grid % fw(s)) * Turb % vis_t_eff(c2) + vis_eff
  end if

  if(c2 < 0) then
    if(Turb % model .eq. LES_SMAGORINSKY    .or.  &
       Turb % model .eq. LES_DYNAMIC        .or.  &
       Turb % model .eq. LES_WALE           .or.  &
       Turb % model .eq. LES_TVM            .or.  &
       Turb % model .eq. HYBRID_LES_RANS    .or.  &
       Turb % model .eq. HYBRID_LES_PRANDTL .or.  &
       Turb % model .eq. K_EPS              .or.  &
       Turb % model .eq. K_EPS_ZETA_F) then
      if(Grid % Bnd_Cond_Type(c2) .eq. WALL .or.  &
         Grid % Bnd_Cond_Type(c2) .eq. WALLFL) then
        vis_eff = Turb % vis_w(c1)
      end if
    end if
  end if

  end subroutine
