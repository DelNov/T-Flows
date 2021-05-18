!==============================================================================!
  subroutine Estimate_Big_And_Small(Grid, big, small)
!------------------------------------------------------------------------------!
!   Estimates "big" and "small" numbers for a given mesh.                      !
!   (They are needed when merging nodes of the mesh.)                          !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type) :: Grid
  real             :: big, small
!-----------------------------------[Locals]-----------------------------------!
  integer :: c, n, l1, l2, n1, n2
  real    :: min_x, min_y, min_z, max_x, max_y, max_z
!==============================================================================!

  !-----------------------------------------------!
  !   Try to work out a reasonable "big" number   !
  !-----------------------------------------------!
  min_x = HUGE; max_x = -HUGE
  min_y = HUGE; max_y = -HUGE
  min_z = HUGE; max_z = -HUGE
  do n = 1, Grid % n_nodes
    min_x = min(min_x, Grid % xn(n));  max_x = max(max_x, Grid % xn(n))
    min_y = min(min_y, Grid % yn(n));  max_y = max(max_y, Grid % yn(n))
    min_z = min(min_z, Grid % zn(n));  max_z = max(max_z, Grid % zn(n))
  end do

  ! Take big to be largest of all dimensions, times roughly 10 (PI**2)
  big = max(max_x-min_x, max_y-min_y, max_z-min_z) * PI**2

  ! One more correction; it did help in some cases
  big = max(big, 1e+6)

  !-------------------------------------------------!
  !   Try to work out a reasonable "small" number   !
  !-------------------------------------------------!
  small = HUGE
  do c = 1, Grid % n_cells
    do l1  = 1, abs(Grid % cells_n_nodes(c))
      n1 = Grid % cells_n(l1, c)
      do l2  = l1+1, abs(Grid % cells_n_nodes(c))
        n2 = Grid % cells_n(l2, c)
        small = min(small, sqrt(   (Grid % xn(n1) - Grid % xn(n2))**2  &
                                 + (Grid % yn(n1) - Grid % yn(n2))**2  &
                                 + (Grid % zn(n1) - Grid % zn(n2))**2 ))
      end do
    end do
  end do

  ! Make small roughly ten times smaller
  small = small / PI**2

!  ! One more correction; it did help in some cases
!  small = min(small, 1e-6)

  end subroutine
