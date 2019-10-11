!==============================================================================!
  subroutine Grad_Mod_Correct_Bad_Cells(grid, phii)
!------------------------------------------------------------------------------!
!   Corrects the pressure gradients in the cells where they cannot             !
!   be computed, the so called "bad" cells.                                    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  real            :: phii(-grid % n_bnd_cells:grid % n_cells)
!-----------------------------------[Locals]-----------------------------------!
  integer :: c, c1, c2, s
!==============================================================================!

  do c = 1, grid % n_cells
    if(bad_cells(c)) then
      phii(c) = 0.0
    end if
  end do 

  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    if(c2 > 0) then
      if(bad_cells(c1)) phii(c1) = phii(c1) + 0.5*phii(c2)
      if(bad_cells(c2)) phii(c2) = phii(c2) + 0.5*phii(c1)
    end if
  end do

  call Comm_Mod_Exchange_Real(grid, phii)

  end subroutine
