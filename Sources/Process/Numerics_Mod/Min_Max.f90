!==============================================================================!
  subroutine Numerics_Mod_Min_Max(phi, phi_min, phi_max)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Var_Type) :: phi
  real           :: phi_min(-phi % pnt_grid % n_bnd_cells:  &
                             phi % pnt_grid % n_cells)
  real           :: phi_max(-phi % pnt_grid % n_bnd_cells:  &
                             phi % pnt_grid % n_cells)
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
