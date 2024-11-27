!==============================================================================!
  subroutine User_Mod_Calculate_Mean(Turb, n0)
!------------------------------------------------------------------------------!
!   User-defined calculation of time-averaged values.                          !
!   Follow Turb_Mod/Calculate_Mean.f90 as guideline for this function.         !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Turb_Type), target :: Turb
  integer, intent(in)     :: n0
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type), pointer :: Flow
  type(Grid_Type),  pointer :: Grid
  integer                   :: c, n
!==============================================================================!

  if(.not. Turb % statistics) return

  n = Time % Curr_Dt() - n0

  ! Take aliases
  Flow => Turb % pnt_flow
  Grid => Flow % pnt_grid

  if(n > -1) then
    do c = -Grid % n_bnd_cells, Grid % n_cells

      !-----------------!
      !   Mean values   !
      !-----------------!

    end do
  end if

  end subroutine
