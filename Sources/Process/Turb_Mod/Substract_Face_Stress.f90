!==============================================================================!
  subroutine Substract_Face_Stress(Turb, ui_si, ui_di, ui_c1, ui_c2,  &
                                         a_fc, b, s)
!------------------------------------------------------------------------------!
!   Computes turbulent stress on a cell face for all turbulence models.        !
!   It is called from Compute_Momentum, while discretizing diffusion terms.    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Turb_Type), target :: Turb
  real                     :: ui_si, ui_di, ui_c1, ui_c2, a_fc
  real                     :: b(Turb % pnt_grid % n_cells)
  integer                  :: s
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: Grid
  real                      :: a0, f_ex, f_im, vis_tur
  integer                   :: c1, c2
!==============================================================================!

  ! Take alias
  Grid => Turb % pnt_grid

  c1 = Grid % faces_c(1,s)
  c2 = Grid % faces_c(2,s)

  if(Turb % model .eq. RSM_MANCEAU_HANJALIC .or.  &
     Turb % model .eq. RSM_HANJALIC_JAKIRLIC) then
    if(Turb % model_variant .ne. STABILIZED) then

      vis_tur =     (Grid % fw(s)  * Turb % vis_t(c1)  &
              + (1.0-Grid % fw(s)) * Turb % vis_t(c2))

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
