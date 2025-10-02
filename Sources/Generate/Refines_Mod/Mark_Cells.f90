!==============================================================================!
  subroutine Refines_Mod_Mark_Cells(ref, Grid)
!------------------------------------------------------------------------------!
!>  Marks specific regions of a domain for local refinement and then refines
!>  the grid accordingly.
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * Initializes all cells as unmarked for refinement.                        !
!   * Iterates through each refinement level and region, checking if each cell !
!     in the grid falls within the defined refinement region.                  !
!   * Uses the shape parameters (ellipsoid, rectangle, or plane) to determine  !
!     whether a cell is within the refinement region.                          !
!   * Calls Refines_Mod_Refine_Marked_Cells to refine cells that are marked    !
!     for refinement.                                                          !
!   * Resets the marked cells after each refinement level is processed.        !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Refines_Type) :: ref   !! type holding information on refinement
  type(Grid_Type)    :: Grid  !! grid being generated (refined here)
!-----------------------------------[Locals]-----------------------------------!
  integer :: c, lev, reg, n1, n2, n3, n4, n5, n6, n7, n8
  real    :: x1, y1, z1, x8, y8, z8, x0, y0, z0
!==============================================================================!

  ! Set no cell for refinement, intially
  ref % cell_marked(:) = .false.

  do lev = 1, ref % n_levels

    do reg = 1, ref % n_ranges(lev)

      x1 = ref % range(lev,reg) % pnt(1) % x
      y1 = ref % range(lev,reg) % pnt(1) % y
      z1 = ref % range(lev,reg) % pnt(1) % z
      x8 = ref % range(lev,reg) % pnt(2) % x
      y8 = ref % range(lev,reg) % pnt(2) % y
      z8 = ref % range(lev,reg) % pnt(2) % z

      do c = 1, Grid % n_cells
        n1 = Grid % cells_n(1,c)
        n2 = Grid % cells_n(2,c)
        n3 = Grid % cells_n(3,c)
        n4 = Grid % cells_n(4,c)
        n5 = Grid % cells_n(5,c)
        n6 = Grid % cells_n(6,c)
        n7 = Grid % cells_n(7,c)
        n8 = Grid % cells_n(8,c)

        x0=1.25e-1*(Grid % xn(n1) + Grid % xn(n2) +   &
                    Grid % xn(n3) + Grid % xn(n4) +   &
                    Grid % xn(n5) + Grid % xn(n6) +   &
                    Grid % xn(n7) + Grid % xn(n8))
        y0=1.25e-1*(Grid % yn(n1) + Grid % yn(n2) +   &
                    Grid % yn(n3) + Grid % yn(n4) +   &
                    Grid % yn(n5) + Grid % yn(n6) +   &
                    Grid % yn(n7) + Grid % yn(n8))
        z0=1.25e-1*(Grid % zn(n1) + Grid % zn(n2) +   &
                    Grid % zn(n3) + Grid % zn(n4) +   &
                    Grid % zn(n5) + Grid % zn(n6) +   &
                    Grid % zn(n7) + Grid % zn(n8))

        if(ref % range(lev,reg) % shape .eq. ELIPSOID) then
          if(  ( ((x1-x0)/x8)**2 +                                  &
                 ((y1-y0)/y8)**2 +                                  &
                 ((z1-z0)/z8)**2)  < 1.0 ) then
            ref % cell_marked(c) = .true.
          end if
        else if(ref % range(lev,reg) % shape .eq. RECTANGLE) then
          if( (x1 < x0) .and. (x0 < x8) .and.                     &
              (y1 < y0) .and. (y0 < y8) .and.                     &
              (z1 < z0) .and. (z0 < z8) ) then
            ref % cell_marked(c) = .true.
          end if
        else if(ref % range(lev,reg) % shape .eq. PLANE) then
          if( (x0-x1)*x8+(y0-y1)*y8+(z0-z1)*z8   >  0.0 ) then
            ref % cell_marked(c) = .true.
          end if
        end if
      end do   ! cells

    end do   ! reg

    call Refines_Mod_Refine_Marked_Cells(ref, Grid, lev)

    ref % cell_marked(:) = .false.

  end do  ! lev

  end subroutine
