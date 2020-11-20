!==============================================================================!
  subroutine Multiphase_Averaging(mult, var)
!------------------------------------------------------------------------------!
!   Averages something in case of phase change.                                !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Multiphase_Type), target :: mult
  type(Var_Type),        target :: var
!-----------------------------------[Locals]-----------------------------------!
  integer :: i, j
  real    :: avg
!==============================================================================!

  if(mult % phase_change) then
    do i = 1, size(mult % avg_cells, 1)
      avg = 0.0
      do j = 1, size(mult % avg_cells, 2)
        avg = avg + var % n(mult % avg_cells(i,j))
      end do

      avg = avg / real(size(mult % avg_cells, 2))

      do j = 1, size(mult % avg_cells, 2)
        var % n(mult % avg_cells(i,j)) = avg
      end do
    end do
  end if

  end subroutine
