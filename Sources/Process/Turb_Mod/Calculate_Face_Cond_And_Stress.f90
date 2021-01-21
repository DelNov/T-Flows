!==============================================================================!
  subroutine Turb_Mod_Calculate_Face_Cond_And_Stress(turb, con_eff, t_stress, s)
!------------------------------------------------------------------------------!
!   Computes turbulent conductivity on a cell face for all turbulence models.  !
!   It is called from Compute_Energy, while discretizing diffusion terms.      !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Turb_Type), target :: turb
  real                    :: con_eff
  real                    :: t_stress
  integer                 :: s
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: grid
  type(Field_Type), pointer :: flow
  type(Var_Type),   pointer :: t
  type(Var_Type),   pointer :: ut, vt, wt
  integer                   :: c1, c2
  real                      :: pr_t1, pr_t2, pr_tf, con_mol, con_turb
  real                      :: cap_dens_c1, cap_dens_c2, tx_f, ty_f, tz_f
  real                      :: ut_cap_dens, vt_cap_dens, wt_cap_dens
!==============================================================================!

  ! Take alias
  flow => turb % pnt_flow
  grid => turb % pnt_grid
  t    => flow % t

  call Turb_Mod_Alias_Heat_Fluxes(turb, ut, vt, wt)

  c1 = grid % faces_c(1,s)
  c2 = grid % faces_c(2,s)

  !------------------------------!
  !   Turbulent Prandtl number   !
  !------------------------------!
  pr_tf = pr_t
  if(turb % model .ne. LES_SMAGORINSKY     .and.  &
     turb % model .ne. LES_DYNAMIC         .and.  &
     turb % model .ne. HYBRID_LES_PRANDTL  .and.  &
     turb % model .ne. LES_WALE            .and.  &
     turb % model .ne. NO_TURBULENCE_MODEL .and.  &
     turb % model .ne. DNS) then
    pr_t1 = Turb_Mod_Prandtl_Number(turb, c1)
    pr_t2 = Turb_Mod_Prandtl_Number(turb, c2)
    pr_tf = grid % fw(s) * pr_t1 + (1.0-grid % fw(s)) * pr_t2
  end if

  !------------------------------------------------------!
  !   Molecular conductivity (without turbulent parts)   !
  !------------------------------------------------------!
  con_mol  =      grid % fw(s)  * flow % conductivity(c1)  &
           + (1.0-grid % fw(s)) * flow % conductivity(c2)
  con_turb = 0.0

  !-----------------------------------------------------------------!
  !   Compute turbulent conductivity for various turbulent models   !
  !-----------------------------------------------------------------!
  if(turb % model .ne. NO_TURBULENCE_MODEL .and.  &
     turb % model .ne. DNS) then
    con_turb  = grid % fw(s) *flow % capacity(c1) * turb % vis_t(c1) / pr_tf  &
       +   (1.0-grid % fw(s))*flow % capacity(c2) * turb % vis_t(c2) / pr_tf
  end if

  if(turb % model .eq. HYBRID_LES_RANS) then
    con_turb  = grid % fw(s)* flow % capacity(c1)*turb % vis_t_eff(c1)/pr_tf  &
         + (1.0-grid % fw(s))*flow % capacity(c2)*turb % vis_t_eff(c2)/pr_tf
  end if

  !-----------------------------------!
  !   Update effective conductivity   !
  !-----------------------------------!
  con_eff = con_mol + con_turb

  ! Effective conductivity at walls
  if(turb % model .eq. K_EPS        .or.  &
     turb % model .eq. K_EPS_ZETA_F .or.  &
     turb % model .eq. HYBRID_LES_RANS) then
    if(c2 < 0) then
      if(Var_Mod_Bnd_Cond_Type(t, c2) .eq. WALL .or.  &
         Var_Mod_Bnd_Cond_Type(t, c2) .eq. WALLFL) then
        con_eff = turb % con_w(c1)
      end if
    end if
  end if

  !---------------------------!
  !   Turbulent heat fluxes   !
  !---------------------------!
  t_stress = 0.0
  if(turb % model .eq. RSM_HANJALIC_JAKIRLIC .or.  &
     turb % model .eq. RSM_MANCEAU_HANJALIC) then

    ! Gradients on the cell face (fw corrects situation close to the wall)
    tx_f = grid % fw(s) * t % x(c1) + (1.0-grid % fw(s)) * t % x(c2)
    ty_f = grid % fw(s) * t % y(c1) + (1.0-grid % fw(s)) * t % y(c2)
    tz_f = grid % fw(s) * t % z(c1) + (1.0-grid % fw(s)) * t % z(c2)

    cap_dens_c1 = flow % capacity(c1) * flow % density(c1)
    cap_dens_c2 = flow % capacity(c2) * flow % density(c2)

    ! Turbulent heat fluxes according to GGDH scheme
    ! (first line is GGDH, second line is SGDH substratced
    ut_cap_dens =  (    grid % fw(s)  * ut % n(c1) * cap_dens_c1    &
                +  (1.0-grid % fw(s)) * ut % n(c2) * cap_dens_c2)
    vt_cap_dens =  (    grid % fw(s)  * vt % n(c1) * cap_dens_c1    &
                +  (1.0-grid % fw(s)) * vt % n(c2) * cap_dens_c2)
    wt_cap_dens =  (    grid % fw(s)  * wt % n(c1) * cap_dens_c1    &
                +  (1.0-grid % fw(s)) * wt % n(c2) * cap_dens_c2)
    t_stress = - (  ut_cap_dens * grid % sx(s)                      &
                  + vt_cap_dens * grid % sy(s)                      &
                  + wt_cap_dens * grid % sz(s) )                    &
                  - (con_turb * (  tx_f * grid % sx(s)              &
                                 + ty_f * grid % sy(s)              &
                                 + tz_f * grid % sz(s)) )

  end if  ! if models are of RSM type

  end subroutine
