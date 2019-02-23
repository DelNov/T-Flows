!==============================================================================!
  subroutine User_Mod_Calculate_Mean(flow, n0, n1)
!------------------------------------------------------------------------------!
!   User-defined calculation of time-averaged values.                          !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type), target :: flow
  integer                  :: n0, n1
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: grid
  integer                  :: c, n
!==============================================================================!

  grid => flow % pnt_grid

  n = n1-n0

  if(n  > -1) then
    do c = -grid % n_bnd_cells, grid % n_cells

      !-----------------!
      !   Mean values   !
      !-----------------!

    end do 
  end if

  end subroutine
