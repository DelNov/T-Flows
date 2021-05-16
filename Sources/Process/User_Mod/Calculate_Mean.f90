!==============================================================================!
  subroutine User_Mod_Calculate_Mean(turb, n0, n1)
!------------------------------------------------------------------------------!
!   User-defined calculation of time-averaged values.                          !
!   Follow Turb_Mod/Calculate_Mean.f90 as guideline for this function.         !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Turb_Type), target  :: turb
  integer, intent(in)      :: n0
  integer, intent(in)      :: n1
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type), pointer :: Flow
  type(Grid_Type),  pointer :: grid
  integer                   :: c, n
!==============================================================================!

  if(.not. turb % statistics) return

  n = n1 - n0

  ! Take aliases
  Flow => turb % pnt_flow
  grid => Flow % pnt_grid

  if(n > -1) then
    do c = -grid % n_bnd_cells, grid % n_cells

      !-----------------!
      !   Mean values   !
      !-----------------!

    end do
  end if

  end subroutine
