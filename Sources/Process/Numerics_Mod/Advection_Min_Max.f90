!==============================================================================!
  subroutine Numerics_Mod_Advection_Min_Max(phi)
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Comm_Mod
  use Grid_Mod, only: Grid_Type
  use Var_Mod,  only: Var_Type
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Var_Type) :: phi
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: grid
  integer                  :: c1, c2, s
!==============================================================================!

  ! Take alias to grid
  grid => phi % pnt_grid

  phi % min(:) = phi % n(:)
  phi % max(:) = phi % n(:)

  do s = 1, grid % n_faces
    c1 = grid % faces_c(1, s)
    c2 = grid % faces_c(2, s)

    phi % min(c1) = min(phi % min(c1), phi % n(c2))
    phi % min(c2) = min(phi % min(c2), phi % n(c1))

    phi % max(c1) = max(phi % max(c1), phi % n(c2))
    phi % max(c2) = max(phi % max(c2), phi % n(c1))

  end do

  call Comm_Mod_Exchange_Real(grid, phi % min)
  call Comm_Mod_Exchange_Real(grid, phi % max)

  end subroutine
