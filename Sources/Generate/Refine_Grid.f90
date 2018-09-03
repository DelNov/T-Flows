!==============================================================================!
  subroutine Refine_Grid(grid)
!------------------------------------------------------------------------------!
!   Mark the region of the domain for local refinement and refine the grid!    !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Gen_Mod
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer :: c, lev, reg, n1, n2, n3, n4, n5, n6, n7, n8
  real    :: x1, y1, z1, x8, y8, z8, x0, y0, z0
!==============================================================================!

  ! Set no cell for refinement, intially
  cell_marked = .false.

  do lev = 1,n_refine_levels

    do reg = 1,n_refined_regions(lev) 

      x1=refined_regions(lev,reg,1)
      y1=refined_regions(lev,reg,2)
      z1=refined_regions(lev,reg,3)
      x8=refined_regions(lev,reg,4)
      y8=refined_regions(lev,reg,5)
      z8=refined_regions(lev,reg,6)

      do c = 1, grid % n_cells
        n1 = grid % cells_n(1,c)
        n2 = grid % cells_n(2,c)
        n3 = grid % cells_n(3,c)
        n4 = grid % cells_n(4,c)
        n5 = grid % cells_n(5,c)
        n6 = grid % cells_n(6,c)
        n7 = grid % cells_n(7,c)
        n8 = grid % cells_n(8,c)

        x0=1.25e-1*(grid % xn(n1) + grid % xn(n2) +   &
                    grid % xn(n3) + grid % xn(n4) +   &
                    grid % xn(n5) + grid % xn(n6) +   &
                    grid % xn(n7) + grid % xn(n8))
        y0=1.25e-1*(grid % yn(n1) + grid % yn(n2) +   &
                    grid % yn(n3) + grid % yn(n4) +   &
                    grid % yn(n5) + grid % yn(n6) +   &
                    grid % yn(n7) + grid % yn(n8))
        z0=1.25e-1*(grid % zn(n1) + grid % zn(n2) +   &
                    grid % zn(n3) + grid % zn(n4) +   &
                    grid % zn(n5) + grid % zn(n6) +   &
                    grid % zn(n7) + grid % zn(n8))

        if(refined_regions(lev,reg,0) .eq. ELIPSOID) then
          if(  ( ((x1-x0)/x8)**2 +                                  &
                 ((y1-y0)/y8)**2 +                                  &
                 ((z1-z0)/z8)**2)  < 1.0 ) then
            cell_marked(c) = .true.
          end if
        else if(refined_regions(lev,reg,0) .eq. RECTANGLE) then 
          if( (x1  < x0) .and. (x0  < x8) .and.                     &
              (y1  < y0) .and. (y0  < y8) .and.                     &
              (z1  < z0) .and. (z0  < z8) ) then
            cell_marked(c) = .true.
          end if
        else if(refined_regions(lev,reg,0) .eq. PLANE) then 
          if( (x0-x1)*x8+(y0-y1)*y8+(z0-z1)*z8   >  0.0 ) then
            cell_marked(c) = .true.
          end if
        end if 
      end do   ! cells

    end do   ! reg

    call Refine_Marked_Cells(grid, lev)

    cell_marked = .false.

  end do  ! lev

  end subroutine
