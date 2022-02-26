!==============================================================================!
  subroutine Face_Stress(Turb, ui, f_stress, s)
!------------------------------------------------------------------------------!
!   Computes turbulent stress on a cell face for all turbulence models.        !
!   It is called from Compute_Momentum, while discretizing diffusion terms.    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Turb_Type), target :: Turb
  real                     :: f_stress
  integer                  :: s
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: Grid
  type(Var_Type),   pointer :: ui
  type(Var_Type),   pointer :: uu, vv, ww, uv, uw, vw
  integer                   :: c1, c2
  real                      :: uu_f, vv_f, ww_f, uv_f, uw_f, vw_f
!==============================================================================!

  ! Take alias
  Grid => Turb % pnt_grid
  call Turb % Alias_Stresses(uu, vv, ww, uv, uw, vw)

  c1 = Grid % faces_c(1,s)
  c2 = Grid % faces_c(2,s)

  ! Add influence of Re stresses
  if(Turb % model .eq. RSM_MANCEAU_HANJALIC .or.  &
     Turb % model .eq. RSM_HANJALIC_JAKIRLIC) then
    if(Turb % model_variant .ne. STABILIZED) then
      if(ui % name .eq. 'U') then
        uu_f = Grid % fw(s) * uu % n(c1) + (1.0-Grid % fw(s)) * uu % n(c2)
        uv_f = Grid % fw(s) * uv % n(c1) + (1.0-Grid % fw(s)) * uv % n(c2)
        uw_f = Grid % fw(s) * uw % n(c1) + (1.0-Grid % fw(s)) * uw % n(c2)
        f_stress = - (  uu_f * Grid % sx(s)  &
                      + uv_f * Grid % sy(s)  &
                      + uw_f * Grid % sz(s) )
      else if(ui % name .eq. 'V') then
        uv_f = Grid % fw(s) * uv % n(c1) + (1.0-Grid % fw(s)) * uv % n(c2)
        vv_f = Grid % fw(s) * vv % n(c1) + (1.0-Grid % fw(s)) * vv % n(c2)
        vw_f = Grid % fw(s) * vw % n(c1) + (1.0-Grid % fw(s)) * vw % n(c2)
        f_stress =  - (  uv_f * Grid % sx(s)  &
                       + vv_f * Grid % sy(s)  &
                       + vw_f * Grid % sz(s) )
      else if(ui % name .eq. 'W') then
        uw_f = Grid % fw(s) * uw % n(c1) + (1.0-Grid % fw(s)) * uw % n(c2)
        vw_f = Grid % fw(s) * vw % n(c1) + (1.0-Grid % fw(s)) * vw % n(c2)
        ww_f = Grid % fw(s) * ww % n(c1) + (1.0-Grid % fw(s)) * ww % n(c2)
        f_stress =  - (  uw_f * Grid % sx(s)  &
                       + vw_f * Grid % sy(s)  &
                       + ww_f * Grid % sz(s) )
      end if
    end if
  end if

  end subroutine
