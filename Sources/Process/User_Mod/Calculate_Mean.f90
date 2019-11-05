!==============================================================================!
  subroutine User_Mod_Calculate_Mean(turb, n0, n1)
!------------------------------------------------------------------------------!
!   User-defined calculation of time-averaged values.                          !
!   Follow Turb_Mod/Calculate_Mean.f90 as guideline for this function.         !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Turb_Type), target  :: turb
  integer                  :: n0, n1
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type), pointer :: flow
  type(Grid_Type),  pointer :: grid
  integer                   :: c, n
!==============================================================================!

  if(.not. turbulence_statistics) return

  n = n1-n0

  ! Take aliases
  flow => turb % pnt_flow
  grid => flow % pnt_grid

  if(n > -1) then
    do c = -grid % n_bnd_cells, grid % n_cells

      !-----------------!
      !   Mean values   !
      !-----------------!

    end do
  end if

  end subroutine
