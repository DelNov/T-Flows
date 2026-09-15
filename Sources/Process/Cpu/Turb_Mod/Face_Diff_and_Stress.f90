!==============================================================================!
  subroutine Face_Diff_And_Stress(Turb, dif_eff, phi_stress, s, sc)
!------------------------------------------------------------------------------!
!   Computes turbulent diffusivity on a cell face for all turbulence models.   !
!   It is called from Compute_Scalar, while discretizing diffusion terms.      !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Turb_Type), target :: Turb
  real,        intent(out) :: dif_eff
  real,        intent(out) :: phi_stress
  integer,     intent(in)  :: s
  integer,     intent(in)  :: sc
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: Grid
  type(Field_Type), pointer :: Flow
  type(Var_Type),   pointer :: phi
  integer                   :: c1, c2
  real                      :: dif_mol, dif_turb
  real                      :: lf1, lf2, alpha_d1, alpha_d2
  real                      :: sc_t1, sc_t2
  real                      :: phix_f, phiy_f, phiz_f
  real                      :: uc_f, vc_f, wc_f
!==============================================================================!

  ! Take alias
  Flow => Turb % pnt_flow
  Grid => Turb % pnt_grid
  phi  => Flow % scalar(sc)

  phi_stress = 0.0

  c1 = Grid % faces_c(1,s)
  c2 = Grid % faces_c(2,s)

  !------------------------------------------------------!
  !   Molecular diffusivity (without turbulent parts)    !
  !------------------------------------------------------!
  dif_mol = Grid % fw(s)  * Flow % diffusivity(c1)  &
     + (1.0-Grid % fw(s)) * Flow % diffusivity(c2)

  dif_turb = 0.0

  !-----------------------------------------------------------------!
  !   Compute turbulent diffusivity for various turbulent models    !
  !-----------------------------------------------------------------!

  if(Turb % model .ne. NO_TURBULENCE_MODEL .and.  &
     Turb % model .ne. DNS) then
    dif_turb = Grid % fw(s)  * Turb % vis_t(c1) / sc_t  &
        + (1.0-Grid % fw(s)) * Turb % vis_t(c2) / sc_t
  end if

  if(Turb % model .eq. HYBRID_LES_RANS) then

    ! Use the same local RANS/LES switching criterion as in
    ! Vis_T_K_Eps_Zeta_F.  RANS cells use the prescribed turbulent
    ! Schmidt number, while LES cells use Sc_t = 0.7.
    lf1      = Grid % vol(c1)**ONE_THIRD
    alpha_d1 = Turb % kappa * Grid % wall_dist(c1) / lf1

    if(alpha_d1 > Turb % c_hyb) then
      sc_t1 = 0.7
    else
      sc_t1 = sc_t
    end if

    if(c2 > 0) then
      lf2      = Grid % vol(c2)**ONE_THIRD
      alpha_d2 = Turb % kappa * Grid % wall_dist(c2) / lf2

      if(alpha_d2 > Turb % c_hyb) then
        sc_t2 = 0.7
      else
        sc_t2 = sc_t
      end if
    else
      sc_t2 = sc_t1
    end if

    ! Turb % vis_t already contains the locally selected RANS/LES viscosity.
    dif_turb = Grid % fw(s) * Turb % vis_t(c1) / sc_t1  &
       + (1.0-Grid % fw(s)) * Turb % vis_t(c2) / sc_t2
  end if

  !-----------------------------------!
  !   Update effective diffusivity    !
  !-----------------------------------!
  dif_eff = dif_mol + dif_turb

  ! Effective diffusivity at walls
  if(Turb % model .eq. K_EPS           .or.  &
     Turb % model .eq. K_EPS_ZETA_F    .or.  &
     Turb % model .eq. SPALART_ALLMARAS.or.  &
     Turb % model .eq. DES_SPALART     .or.  &
     Turb % model .eq. LES_SMAGORINSKY .or.  &
     Turb % model .eq. LES_DYNAMIC     .or.  &
     Turb % model .eq. LES_WALE        .or.  &
     Turb % model .eq. LES_TVM         .or.  &
     Turb % model .eq. HYBRID_LES_RANS) then
    if(c2 < 0) then
      if(Var_Mod_Bnd_Cond_Type(phi, c2) .eq. WALL .or.  &
         Var_Mod_Bnd_Cond_Type(phi, c2) .eq. WALLFL) then
        dif_eff = Turb % diff_w(c1)
      end if
    end if
  end if

  !---------------------------!
  !   Turbulent scalar fluxes !
  !---------------------------!

  if(Turb % scalar_flux_model .eq. GGDH .or. &
     Turb % scalar_flux_model .eq. AFM) then

    ! Gradients on the cell face (fw corrects situation close to the wall)
    phix_f = Grid % fw(s) * phi % x(c1) + (1.0-Grid % fw(s)) * phi % x(c2)
    phiy_f = Grid % fw(s) * phi % y(c1) + (1.0-Grid % fw(s)) * phi % y(c2)
    phiz_f = Grid % fw(s) * phi % z(c1) + (1.0-Grid % fw(s)) * phi % z(c2)

    ! Turbulent heat fluxes according to GGDH scheme
    ! (first line is GGDH, second line is SGDH substratced
    uc_f =  (    Grid % fw(s)  * Turb % uc(c1) * Flow % density(c1)   &
         +  (1.0-Grid % fw(s)) * Turb % uc(c2) * Flow % density(c2))
    vc_f =  (    Grid % fw(s)  * Turb % vc(c1) * Flow % density(c1)   &
         +  (1.0-Grid % fw(s)) * Turb % vc(c2) * Flow % density(c2))
    wc_f =  (    Grid % fw(s)  * Turb % wc(c1) * Flow % density(c1)   &
         +  (1.0-Grid % fw(s)) * Turb % wc(c2) * Flow % density(c2))

    phi_stress = - (  uc_f * Grid % sx(s)                      &
                    + vc_f * Grid % sy(s)                      &
                    + wc_f * Grid % sz(s) )                    &
                  - (dif_turb * (  phix_f * Grid % sx(s)              &
                                 + phiy_f * Grid % sy(s)              &
                                 + phiz_f * Grid % sz(s)) )

    if(Grid % cell_near_wall(c1).or.Grid % cell_near_wall(c2)) then
      if( Turb % y_plus(c1) > 11.0 .or. Turb % y_plus(c2) > 11.0 ) then

        phi_stress = 0.0

      end if
    end if

  end if  ! if model is AFM of GGDH

  end subroutine
