!==============================================================================!
  subroutine Turb_Mod_Face_Diff_And_Stress(turb, dif_eff, phi_stress, s, sc)
!------------------------------------------------------------------------------!
!   Computes turbulent diffusivity on a cell face for all turbulence models.   !
!   It is called from Compute_Scalar, while discretizing diffusion terms.      !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Turb_Type), target :: turb
  real                    :: dif_eff
  real                    :: phi_stress
  integer                 :: s
  integer, intent(in)     :: sc
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: Grid
  type(Field_Type), pointer :: Flow
  type(Var_Type),   pointer :: phi
  integer                   :: c1, c2
  real                      :: dif_mol, dif_turb
  real                      :: phix_f, phiy_f, phiz_f
  real                      :: uc_f, vc_f, wc_f
!==============================================================================!

  ! Take alias
  Flow => turb % pnt_flow
  Grid => turb % pnt_grid
  phi  => Flow % scalar(sc)

  c1 = Grid % faces_c(1,s)
  c2 = Grid % faces_c(2,s)

  !------------------------------------------------------!
  !   Molecular diffusivity (without turbulent parts)    !
  !------------------------------------------------------!
  dif_mol  = Flow % diffusivity 
             
  dif_turb = 0.0

  !-----------------------------------------------------------------!
  !   Compute turbulent diffusivity for various turbulent models    !
  !-----------------------------------------------------------------!

  if(turb % model .ne. NO_TURBULENCE_MODEL .and.  &
     turb % model .ne. DNS) then
    dif_turb = Grid % fw(s) * turb % vis_t(c1) / sc_t  &
       + (1.0-Grid % fw(s)) * turb % vis_t(c2) / sc_t
  end if

  if(turb % model .eq. HYBRID_LES_RANS) then
    dif_turb = Grid % fw(s) * turb % vis_t_eff(c1) / sc_t  &
       + (1.0-Grid % fw(s)) * turb % vis_t_eff(c2) / sc_t
  end if

  !-----------------------------------!
  !   Update effective diffusivity    !
  !-----------------------------------!
  dif_eff = dif_mol + dif_turb 

  ! Effective diffusivity at walls
  if(turb % model .eq. K_EPS        .or.  &
     turb % model .eq. K_EPS_ZETA_F .or.  &
     turb % model .eq. HYBRID_LES_RANS) then
    if(c2 < 0) then
      if(Var_Mod_Bnd_Cond_Type(phi, c2) .eq. WALL .or.  &
         Var_Mod_Bnd_Cond_Type(phi, c2) .eq. WALLFL) then
        dif_eff = turb % diff_w(c1)
      end if
    end if
  end if

  !---------------------------!
  !   Turbulent scalar fluxes !
  !---------------------------!

  if(turb % scalar_flux_model .eq. GGDH .or. &
     turb % scalar_flux_model .eq. AFM) then 

    ! Gradients on the cell face (fw corrects situation close to the wall)
    phix_f = Grid % fw(s) * phi % x(c1) + (1.0-Grid % fw(s)) * phi % x(c2)
    phiy_f = Grid % fw(s) * phi % y(c1) + (1.0-Grid % fw(s)) * phi % y(c2)
    phiz_f = Grid % fw(s) * phi % z(c1) + (1.0-Grid % fw(s)) * phi % z(c2)

    ! Turbulent heat fluxes according to GGDH scheme
    ! (first line is GGDH, second line is SGDH substratced
    uc_f =  (    Grid % fw(s)  * turb % uc(c1) * Flow % density(c1)   &
         +  (1.0-Grid % fw(s)) * turb % uc(c2) * Flow % density(c2))
    vc_f =  (    Grid % fw(s)  * turb % vc(c1) * Flow % density(c1)   &
         +  (1.0-Grid % fw(s)) * turb % vc(c2) * Flow % density(c2))
    wc_f =  (    Grid % fw(s)  * turb % wc(c1) * Flow % density(c1)   &
         +  (1.0-Grid % fw(s)) * turb % wc(c2) * Flow % density(c2))

    phi_stress = - (  uc_f * Grid % sx(s)                      &
                    + vc_f * Grid % sy(s)                      &
                    + wc_f * Grid % sz(s) )                    &
                  - (dif_turb * (  phix_f * Grid % sx(s)              &
                                 + phiy_f * Grid % sy(s)              &
                                 + phiz_f * Grid % sz(s)) )
  
    if(Grid % cell_near_wall(c1).or.Grid % cell_near_wall(c2)) then
      if( turb % y_plus(c1) > 11.0 .or. turb % y_plus(c2) > 11.0 ) then
 
        phi_stress = 0.0
  
      end if
    end if
    
  end if  ! if model is AFM of GGDH 

  end subroutine
