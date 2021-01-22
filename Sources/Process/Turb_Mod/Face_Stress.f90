!==============================================================================!
  subroutine Turb_Mod_Face_Stress(turb, ui, f_stress, s)
!------------------------------------------------------------------------------!
!   Computes turbulent stress on a cell face for all turbulence models.        !
!   It is called from Compute_Momentum, while discretizing diffusion terms.    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Turb_Type), target :: turb
  real                    :: f_stress
  integer                 :: s
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: grid
  type(Var_Type),   pointer :: ui
  type(Var_Type),   pointer :: uu, vv, ww, uv, uw, vw
  integer                   :: c1, c2
  real                      :: uu_f, vv_f, ww_f, uv_f, uw_f, vw_f
!==============================================================================!

  ! Take alias
  grid => turb % pnt_grid
  call Turb_Mod_Alias_Stresses(turb, uu, vv, ww, uv, uw, vw)

  c1 = grid % faces_c(1,s)
  c2 = grid % faces_c(2,s)

  ! Add influence of Re stresses
  if(turb % model .eq. RSM_MANCEAU_HANJALIC .or.  &
     turb % model .eq. RSM_HANJALIC_JAKIRLIC) then
    if(turb % model_variant .ne. STABILIZED) then
      if(ui % name .eq. 'U') then
        uu_f = grid % fw(s) * uu % n(c1) + (1.0-grid % fw(s)) * uu % n(c2)
        uv_f = grid % fw(s) * uv % n(c1) + (1.0-grid % fw(s)) * uv % n(c2)
        uw_f = grid % fw(s) * uw % n(c1) + (1.0-grid % fw(s)) * uw % n(c2)
        f_stress = - (  uu_f * grid % sx(s)  &
                      + uv_f * grid % sy(s)  &
                      + uw_f * grid % sz(s) )
      else if(ui % name .eq. 'V') then
        uv_f = grid % fw(s) * uv % n(c1) + (1.0-grid % fw(s)) * uv % n(c2)
        vv_f = grid % fw(s) * vv % n(c1) + (1.0-grid % fw(s)) * vv % n(c2)
        vw_f = grid % fw(s) * vw % n(c1) + (1.0-grid % fw(s)) * vw % n(c2)
        f_stress =  - (  uv_f * grid % sx(s)  &
                       + vv_f * grid % sy(s)  &
                       + vw_f * grid % sz(s) )
      else if(ui % name .eq. 'W') then
        uw_f = grid % fw(s) * uw % n(c1) + (1.0-grid % fw(s)) * uw % n(c2)
        vw_f = grid % fw(s) * vw % n(c1) + (1.0-grid % fw(s)) * vw % n(c2)
        ww_f = grid % fw(s) * ww % n(c1) + (1.0-grid % fw(s)) * ww % n(c2)
        f_stress =  - (  uw_f * grid % sx(s)  &
                       + vw_f * grid % sy(s)  &
                       + ww_f * grid % sz(s) )
      end if
    end if
  end if

  end subroutine
