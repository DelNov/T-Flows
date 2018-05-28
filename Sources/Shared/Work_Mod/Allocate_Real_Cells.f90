!==============================================================================!
  subroutine Work_Mod_Allocate_Real_Cells(grid, n)
!------------------------------------------------------------------------------!
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type), target :: grid
  integer                 :: n    ! number of arrays
!-----------------------------------[Locals]-----------------------------------!
  integer :: nc, nb
!==============================================================================!

  ! Get number of cells and boundary cells
  nc = grid % n_cells
  nb = grid % n_bnd_cells

  ! Store the pointer to the grid
  pnt_grid => grid

  ! Allocate requested memory
  allocate(r_cell_01(-nb:nc));   r_cell_01 = 0.0;   if(n .eq.  1) return
  allocate(r_cell_02(-nb:nc));   r_cell_02 = 0.0;   if(n .eq.  2) return
  allocate(r_cell_03(-nb:nc));   r_cell_03 = 0.0;   if(n .eq.  3) return
  allocate(r_cell_04(-nb:nc));   r_cell_04 = 0.0;   if(n .eq.  4) return
  allocate(r_cell_05(-nb:nc));   r_cell_05 = 0.0;   if(n .eq.  5) return
  allocate(r_cell_06(-nb:nc));   r_cell_06 = 0.0;   if(n .eq.  6) return
  allocate(r_cell_07(-nb:nc));   r_cell_07 = 0.0;   if(n .eq.  7) return
  allocate(r_cell_08(-nb:nc));   r_cell_08 = 0.0;   if(n .eq.  8) return
  allocate(r_cell_09(-nb:nc));   r_cell_09 = 0.0;   if(n .eq.  9) return
  allocate(r_cell_10(-nb:nc));   r_cell_10 = 0.0;   if(n .eq. 10) return
  allocate(r_cell_11(-nb:nc));   r_cell_11 = 0.0;   if(n .eq. 11) return
  allocate(r_cell_12(-nb:nc));   r_cell_12 = 0.0;   if(n .eq. 12) return
  allocate(r_cell_13(-nb:nc));   r_cell_13 = 0.0;   if(n .eq. 13) return
  allocate(r_cell_14(-nb:nc));   r_cell_14 = 0.0;   if(n .eq. 14) return
  allocate(r_cell_15(-nb:nc));   r_cell_15 = 0.0;   if(n .eq. 15) return
  allocate(r_cell_16(-nb:nc));   r_cell_16 = 0.0;   if(n .eq. 16) return
  allocate(r_cell_17(-nb:nc));   r_cell_17 = 0.0;   if(n .eq. 17) return
  allocate(r_cell_18(-nb:nc));   r_cell_18 = 0.0;   if(n .eq. 18) return
  allocate(r_cell_19(-nb:nc));   r_cell_19 = 0.0;   if(n .eq. 19) return
  allocate(r_cell_20(-nb:nc));   r_cell_20 = 0.0;   if(n .eq. 20) return
  allocate(r_cell_21(-nb:nc));   r_cell_21 = 0.0;   if(n .eq. 21) return
  allocate(r_cell_22(-nb:nc));   r_cell_22 = 0.0;   if(n .eq. 22) return
  allocate(r_cell_23(-nb:nc));   r_cell_23 = 0.0;   if(n .eq. 23) return
  allocate(r_cell_24(-nb:nc));   r_cell_24 = 0.0;   if(n .eq. 24) return
  allocate(r_cell_25(-nb:nc));   r_cell_25 = 0.0;   if(n .eq. 25) return
  allocate(r_cell_26(-nb:nc));   r_cell_26 = 0.0;   if(n .eq. 26) return
  allocate(r_cell_27(-nb:nc));   r_cell_27 = 0.0;   if(n .eq. 27) return
  allocate(r_cell_28(-nb:nc));   r_cell_28 = 0.0;   if(n .eq. 28) return
  allocate(r_cell_29(-nb:nc));   r_cell_29 = 0.0;   if(n .eq. 29) return
  allocate(r_cell_30(-nb:nc));   r_cell_30 = 0.0;   if(n .eq. 30) return

  end subroutine
