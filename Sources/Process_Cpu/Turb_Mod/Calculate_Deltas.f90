!==============================================================================!
  subroutine Calculate_Deltas(Turb)
!------------------------------------------------------------------------------!
!   Calculates h_max, h_min and h_w needed for Spalart Allmaras models.        !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Turb_Type), target :: Turb
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),   pointer :: Grid
  type(Field_Type),  pointer :: Flow
  integer                    :: s, c, c1, c2, nc, nb
  real                       :: d, d1, d2, dx1, dx2, dy1, dy2, dz1, dz2
  real, contiguous,  pointer :: h_w_x(:), h_w_y(:), h_w_z(:)
!==============================================================================!

  if( Turb % model .ne. SPALART_ALLMARAS) return

  call Work % Connect_Real_Cell(h_w_x, h_w_y, h_w_z)

  ! Take aliases
  Flow => Turb % pnt_flow
  Grid => Turb % pnt_grid
  nc = Grid % n_cells
  nb = Grid % n_bnd_cells

  ! Compute gradients of wall distance
  call Flow % Grad(Grid, Grid % wall_dist, h_w_x(-nb:nc),  &
                                           h_w_y(-nb:nc),  &
                                           h_w_z(-nb:nc))

  ! Normalize gradients
  do c = Cells_In_Domain_And_Buffers()
    d = sqrt(  h_w_x(c)**2 + h_w_y(c)**2 + h_w_z(c)**2)
    h_w_x(c) = h_w_x(c) / d
    h_w_y(c) = h_w_y(c) / d
    h_w_z(c) = h_w_z(c) / d
  end do

  ! Initialize all "deltas"
  Turb % h_min(:) = +HUGE
  Turb % h_max(:) = -HUGE
  Turb % h_w  (:) = 0.0

  ! Browse through faces
  do s = 1, Grid % n_faces

    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)

    ! Cell deltas
    d = sqrt(Grid % dx(s)**2 + Grid % dy(s)**2 + Grid % dz(s)**2)

    d1 = d *      Grid % fw(s)
    d2 = d * (1.0-Grid % fw(s))

    Turb % h_max(c1) = max(Turb % h_max(c1), d1)
    Turb % h_max(c2) = max(Turb % h_max(c2), d2)
    Turb % h_min(c1) = min(Turb % h_min(c1), d1)
    Turb % h_min(c2) = min(Turb % h_min(c2), d2)

    ! Wall normal distance
    dx1 = Grid % dx(s) *      Grid % fw(s)
    dy1 = Grid % dy(s) *      Grid % fw(s)
    dz1 = Grid % dz(s) *      Grid % fw(s)
    dx2 = Grid % dx(s) * (1.0-Grid % fw(s))
    dy2 = Grid % dy(s) * (1.0-Grid % fw(s))
    dz2 = Grid % dz(s) * (1.0-Grid % fw(s))

    Turb % h_w(c1) = (Turb % h_w(c1)+ abs(  dx1 * h_w_x(c1)  &
                                          + dy1 * h_w_y(c1)  &
                                          + dz1 * h_w_z(c1)))
    if(c2 > 0) then
      Turb % h_w(c2) = (Turb % h_w(c2)+ abs(  dx2 * h_w_x(c2)  &
                                            + dy2 * h_w_y(c2)  &
                                            + dz2 * h_w_z(c2)))
    end if
  end do

  ! Correct h_max and h_min
  do c = Cells_In_Domain_And_Buffers()
    Turb % h_max(c) = Turb % h_max(c) * 2.0
    Turb % h_min(c) = Turb % h_min(c) * 2.0
  end do

  call Grid % Exchange_Cells_Real(Turb % h_max)
  call Grid % Exchange_Cells_Real(Turb % h_min)
  call Grid % Exchange_Cells_Real(Turb % h_w)

  call Work % Disconnect_Real_Cell(h_w_x, h_w_y, h_w_z)

  end subroutine
