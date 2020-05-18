!==============================================================================!
  subroutine Turb_Mod_Calculate_Deltas(turb)
!------------------------------------------------------------------------------!
!   Calculates h_max, h_min and h_w needed for Spalart Allmaras models.        !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Work_Mod, only: h_w_x => r_cell_01,  &
                      h_w_y => r_cell_02,  &
                      h_w_z => r_cell_03
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Turb_Type), target :: turb
!----------------------------------[Calling]-----------------------------------!
  real Distance
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: grid
  type(Field_Type), pointer :: flow
  integer                   :: s, c, c1, c2, nc, nb
  real                      :: d, d1, d2, dx1, dx2, dy1, dy2, dz1, dz2
!==============================================================================!

  ! Take aliases
  flow => turb % pnt_flow
  grid => turb % pnt_grid
  nc = grid % n_cells
  nb = grid % n_bnd_cells

  ! Compute gradients of wall distance
  call Field_Mod_Grad_Component(flow, grid % wall_dist, 1, h_w_x(-nb:nc))
  call Field_Mod_Grad_Component(flow, grid % wall_dist, 2, h_w_y(-nb:nc))
  call Field_Mod_Grad_Component(flow, grid % wall_dist, 3, h_w_z(-nb:nc))

  ! Normalize gradients
  do c = 1, grid % n_cells
    d = sqrt(  h_w_x(c)**2 + h_w_y(c)**2 + h_w_z(c)**2)
    h_w_x(c) = h_w_x(c) / d
    h_w_y(c) = h_w_y(c) / d
    h_w_z(c) = h_w_z(c) / d
  end do

  ! Initialize all "deltas"
  turb % h_min(:) = +HUGE
  turb % h_max(:) = -HUGE
  turb % h_w  (:) = 0.0

  ! Browse through faces
  do s = 1, grid % n_faces

    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    ! Cell deltas
    d = sqrt(grid % dx(s)**2 + grid % dy(s)**2 + grid % dz(s)**2)

    d1 = d *      grid % fw(s)
    d2 = d * (1.0-grid % fw(s))

    turb % h_max(c1) = max(turb % h_max(c1), d1)
    turb % h_max(c2) = max(turb % h_max(c2), d2)
    turb % h_min(c1) = min(turb % h_min(c1), d1)
    turb % h_min(c2) = min(turb % h_min(c2), d2)

    ! Wall normal distance
    dx1 = grid % dx(s) *      grid % fw(s)
    dy1 = grid % dy(s) *      grid % fw(s)
    dz1 = grid % dz(s) *      grid % fw(s)
    dx2 = grid % dx(s) * (1.0-grid % fw(s))
    dy2 = grid % dy(s) * (1.0-grid % fw(s))
    dz2 = grid % dz(s) * (1.0-grid % fw(s))

    turb % h_w(c1) = (turb % h_w(c1)+ abs(  dx1 * h_w_x(c1)  &
                                          + dy1 * h_w_y(c1)  &
                                          + dz1 * h_w_z(c1)))
    if(c2 > 0) then
      turb % h_w(c2) = (turb % h_w(c2)+ abs(  dx2 * h_w_x(c2)  &
                                            + dy2 * h_w_y(c2)  &
                                            + dz2 * h_w_z(c2)))
    end if
  end do

  ! Correct h_max and h_min
  do c = 1, grid % n_cells
    turb % h_max(c) = turb % h_max(c) * 2.0
    turb % h_min(c) = turb % h_min(c) * 2.0
  end do

  call Grid_Mod_Exchange_Real(grid, turb % h_max)
  call Grid_Mod_Exchange_Real(grid, turb % h_min)
  call Grid_Mod_Exchange_Real(grid, turb % h_w)

  end subroutine
