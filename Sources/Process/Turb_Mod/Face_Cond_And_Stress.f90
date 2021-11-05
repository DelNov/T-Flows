!==============================================================================!
  subroutine Turb_Mod_Face_Cond_And_Stress(turb, con_eff, t_stress, s)
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
  type(Grid_Type),  pointer :: Grid
  type(Field_Type), pointer :: Flow
  type(Var_Type),   pointer :: t, t2
  type(Var_Type),   pointer :: ut, vt, wt
  integer                   :: c1, c2
  real                      :: pr_t1, pr_t2, pr_tf, con_mol, con_turb
  real                      :: cap_dens_c1, cap_dens_c2, tx_f, ty_f, tz_f
  real                      :: ut_cap_dens, vt_cap_dens, wt_cap_dens
!==============================================================================!

  ! Take alias
  Flow => turb % pnt_flow
  Grid => turb % pnt_grid
  t    => Flow % t

  call Turb_Mod_Alias_Heat_Fluxes(turb, ut, vt, wt)
  call Turb_Mod_Alias_T2         (turb, t2)

  c1 = Grid % faces_c(1,s)
  c2 = Grid % faces_c(2,s)

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
    pr_tf = Grid % fw(s) * pr_t1 + (1.0-Grid % fw(s)) * pr_t2
  end if

  !------------------------------------------------------!
  !   Molecular conductivity (without turbulent parts)   !
  !------------------------------------------------------!
  con_mol = 2.0 / (    1.0 / Flow % conductivity(c1)  &
                     + 1.0 / Flow % conductivity(c2) )
  con_turb = 0.0

  !-----------------------------------------------------------------!
  !   Compute turbulent conductivity for various turbulent models   !
  !-----------------------------------------------------------------!
  if(turb % model .ne. NO_TURBULENCE_MODEL .and.  &
     turb % model .ne. DNS) then
    con_turb  = Grid % fw(s) *Flow % capacity(c1) * turb % vis_t(c1) / pr_tf  &
       +   (1.0-Grid % fw(s))*Flow % capacity(c2) * turb % vis_t(c2) / pr_tf
  end if

  if(turb % model .eq. HYBRID_LES_RANS) then
    con_turb  = Grid % fw(s)* Flow % capacity(c1)*turb % vis_t_eff(c1)/pr_tf  &
         + (1.0-Grid % fw(s))*Flow % capacity(c2)*turb % vis_t_eff(c2)/pr_tf
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

  if(turb % heat_flux_model .eq. GGDH .or. &
     turb % heat_flux_model .eq. AFM) then 

    ! Gradients on the cell face (fw corrects situation close to the wall)
    tx_f = Grid % fw(s) * t % x(c1) + (1.0-Grid % fw(s)) * t % x(c2)
    ty_f = Grid % fw(s) * t % y(c1) + (1.0-Grid % fw(s)) * t % y(c2)
    tz_f = Grid % fw(s) * t % z(c1) + (1.0-Grid % fw(s)) * t % z(c2)

    cap_dens_c1 = Flow % capacity(c1) * Flow % density(c1)
    cap_dens_c2 = Flow % capacity(c2) * Flow % density(c2)

    ! Turbulent heat fluxes according to GGDH scheme
    ! (first line is GGDH, second line is SGDH substratced
    ut_cap_dens =  (    Grid % fw(s)  * ut % n(c1) * cap_dens_c1    &
                +  (1.0-Grid % fw(s)) * ut % n(c2) * cap_dens_c2)
    vt_cap_dens =  (    Grid % fw(s)  * vt % n(c1) * cap_dens_c1    &
                +  (1.0-Grid % fw(s)) * vt % n(c2) * cap_dens_c2)
    wt_cap_dens =  (    Grid % fw(s)  * wt % n(c1) * cap_dens_c1    &
                +  (1.0-Grid % fw(s)) * wt % n(c2) * cap_dens_c2)

    t_stress = - (  ut_cap_dens * Grid % sx(s)                      &
                  + vt_cap_dens * Grid % sy(s)                      &
                  + wt_cap_dens * Grid % sz(s) )                    &
                  - (con_turb * (  tx_f * Grid % sx(s)              &
                                 + ty_f * Grid % sy(s)              &
                                 + tz_f * Grid % sz(s)) )

    if(Grid % cell_near_wall(c1).or.Grid % cell_near_wall(c2)) then
      if( turb % y_plus(c1) > 11.0 .or. turb % y_plus(c2) > 11.0 ) then
 
        t_stress = 0.0
  
      end if
    end if
  end if  ! if models are of RSM type

  end subroutine
