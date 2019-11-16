!==============================================================================!
  subroutine Turb_Mod_Substract_Face_Stress(turb, ui_si, ui_di, ui_c1, ui_c2,  &
                                                  a_fc, b, s)
!------------------------------------------------------------------------------!
!   Computes turbulent stress on a cell face for all turbulence models.        !
!   It is called from Compute_Momentum, while discretizing diffusion terms.    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Turb_Type), target :: turb
  real                    :: ui_si, ui_di, ui_c1, ui_c2, a_fc
  real                    :: b(turb % pnt_grid % n_cells)
  integer                 :: s
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: grid
  real                      :: a0, f_ex, f_im, vis_tur
  integer                   :: c1, c2
!==============================================================================!

  ! Take alias
  grid => turb % pnt_grid
! call Turb_Mod_Alias_Stresses(turb, uu, vv, ww, uv, uw, vw)

  c1 = grid % faces_c(1,s)
  c2 = grid % faces_c(2,s)

  if(turbulence_model .eq. RSM_MANCEAU_HANJALIC .or.  &
     turbulence_model .eq. RSM_HANJALIC_JAKIRLIC) then
    if(turbulence_model_variant .ne. STABILIZED) then

      vis_tur =     (grid % fw(s)  * turb % vis_t(c1)  &
              + (1.0-grid % fw(s)) * turb % vis_t(c2))

      f_ex = vis_tur * ui_si

      a0 = a_fc * vis_tur
      f_im = a0 * ui_di

      b(c1) = b(c1) - a0 * (ui_c2 - ui_c1) - f_ex + f_im
      if(c2  > 0) then
        b(c2) = b(c2) + a0 * (ui_c2 - ui_c1) + f_ex - f_im
      end if
    end if
  end if

  end subroutine
