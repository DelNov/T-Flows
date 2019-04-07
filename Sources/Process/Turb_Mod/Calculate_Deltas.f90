!==============================================================================!
  subroutine Calculate_Deltas(grid)
!------------------------------------------------------------------------------!
!   Calculates h_max, needed for Spalart Allmaras and related models.          !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Comm_Mod
  use Grid_Mod
  use Grad_Mod
  use Turb_Mod
  use Work_Mod, only: h_w_x => r_cell_01,  &
                      h_w_y => r_cell_02,  &
                      h_w_z => r_cell_03
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!----------------------------------[Calling]-----------------------------------!
  real Distance
!-----------------------------------[Locals]-----------------------------------!
  integer :: s, c, c1, c2
  real    :: d, d1, d2, dx1, dx2, dy1, dy2, dz1, dz2,  &
                        xc1, xc2, yc1, yc2, zc1, zc2
!==============================================================================!

  ! Compute gradients of wall distance
  call Grad_Mod_Array(grid, grid % wall_dist, h_w_x, h_w_y, h_w_z, .true.)

  ! Normalize gradients
  do c = 1, grid % n_cells
    d = sqrt(  h_w_x(c)**2 + h_w_y(c)**2 + h_w_z(c)**2)
    h_w_x(c) = h_w_x(c) / d
    h_w_y(c) = h_w_y(c) / d
    h_w_z(c) = h_w_z(c) / d
  end do

  ! Initialize all "deltas"
  h_min(:) = +HUGE
  h_max(:) = -HUGE
  h_w  (:) = 0.0

  ! Browse through faces
  do s = 1, grid % n_faces

    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    ! Cell deltas
    d = sqrt(grid % dx(s)**2 + grid % dy(s)**2 + grid % dz(s)**2)

    d1 = d *      grid % fw(s)
    d2 = d * (1.0-grid % fw(s))

    h_max(c1) = max(h_max(c1), d1)
    h_max(c2) = max(h_max(c2), d2)
    h_min(c1) = min(h_min(c1), d1)
    h_min(c2) = min(h_min(c2), d2)

    ! Wall normal distance
    dx1 = grid % dx(s) *      grid % fw(s)
    dy1 = grid % dy(s) *      grid % fw(s)
    dz1 = grid % dz(s) *      grid % fw(s)
    dx2 = grid % dx(s) * (1.0-grid % fw(s))
    dy2 = grid % dy(s) * (1.0-grid % fw(s))
    dz2 = grid % dz(s) * (1.0-grid % fw(s))

    h_w(c1) =    (h_w(c1)+ abs(  dx1 * h_w_x(c1)  &
                               + dy1 * h_w_y(c1)  &
                               + dz1 * h_w_z(c1)))
    if(c2 > 0) then
      h_w(c2) =    (h_w(c2)+ abs(  dx2 * h_w_x(c2)  &
                                 + dy2 * h_w_y(c2)  &
                                 + dz2 * h_w_z(c2)))
    end if
  end do

  ! Correct h_max and h_min
  do c = 1, grid % n_cells
    h_max(c) = h_max(c) * 2.0
    h_min(c) = h_min(c) * 2.0
  end do

  call Comm_Mod_Exchange_Real(grid, h_max)
  call Comm_Mod_Exchange_Real(grid, h_min)
  call Comm_Mod_Exchange_Real(grid, h_w)

  end subroutine
