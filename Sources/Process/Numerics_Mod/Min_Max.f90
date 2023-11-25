!==============================================================================!
  subroutine Numerics_Mod_Min_Max(phi, phi_min, phi_max)
!------------------------------------------------------------------------------!
!>  The Numerics_Mod_Min_Max subroutine in the Numerics_Mod module computes 
!>  the minimum and maximum values of a transported variable among neighbours
!>  of each computational cell.  These bounded values are subsequently used
!>  in Numerics_Mod_Advection_Scheme.
!------------------------------------------------------------------------------!
!   Functionality:                                                             !
!                                                                              !
!   * Extrema calculation: Determines the minimum and maximum values (extrema) !
!     of a variable across the neighboring cells, which is critical for        !
!     interpolating face values with various advection schemes.                !
!   * Face-based comparison: Compares values of the variable at the faces      !
!     between adjacent cells, ensuring that extrema are correctly evaluated.   !
!   * Grid-wide Application: Applies the calculations across all faces of the  !
!     grid, ensuring a comprehensive assessment of the variable's range.       !
!   * Exchange mechanism: Utilizes grid exchange functions to synchronize the  !
!     extrema values across different parts of the grid for parallel runs.     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Var_Type), intent(in)  :: phi  !! transported variable
  real,           intent(out) :: phi_min(-phi % pnt_grid % n_bnd_cells:  &
                                          phi % pnt_grid % n_cells)
                                      !! lower bound array for phi
  real,           intent(out) :: phi_max(-phi % pnt_grid % n_bnd_cells:  &
                                          phi % pnt_grid % n_cells)
                                      !! upper bound array for phi
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
  integer                  :: c1, c2, s
!==============================================================================!

  ! Take alias to Grid
  Grid => phi % pnt_grid

  phi_min(:) = phi % n(:)
  phi_max(:) = phi % n(:)

  do s = 1, Grid % n_faces
    c1 = Grid % faces_c(1, s)
    c2 = Grid % faces_c(2, s)

    phi_min(c1) = min(phi_min(c1), phi % n(c2))
    phi_min(c2) = min(phi_min(c2), phi % n(c1))

    phi_max(c1) = max(phi_max(c1), phi % n(c2))
    phi_max(c2) = max(phi_max(c2), phi % n(c1))
  end do

  call Grid % Exchange_Cells_Real(phi_min)
  call Grid % Exchange_Cells_Real(phi_max)

  end subroutine
