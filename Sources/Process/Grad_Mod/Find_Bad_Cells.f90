!==============================================================================!
  subroutine Grad_Mod_Find_Bad_Cells(grid)
!------------------------------------------------------------------------------!
!   Searches for cells which are "bad" for calculation of pressure gradients.  !
!                                                                              !
!   Practically, these are the tetrahedronal cells with two faces on the       !
!   boundary and two in the domain.                                            ! 
!------------------------------------------------------------------------------!
  use Work_Mod, only: n_good => i_cell_01
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer :: s, c, c1, c2, n_bad
!==============================================================================!

  bad_cells = .false. 
  n_good    = 0

  n_bad = 0

  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    if(c2 < 0) then
      n_good(c1) = n_good(c1) + 1 
    end if  
  end do

  do c = 1, grid % n_cells
    if(n_good(c) .eq. 2) then
      bad_cells(c) = .true.
      n_bad = n_bad + 1
    end if
  end do 

  n_good = 0
  
  call Comm_Mod_Global_Sum_Int(n_bad)

  if(this_proc < 2) print *, '# There are ', n_bad, ' bad cells for gradients.'

  end subroutine
