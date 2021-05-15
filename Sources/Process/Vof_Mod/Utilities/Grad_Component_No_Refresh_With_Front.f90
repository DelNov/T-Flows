!==============================================================================!
  subroutine Grad_Component_No_Refresh_With_Front(Vof, phi, i, phii, phif)
!------------------------------------------------------------------------------!
!   Calculates gradient of generic variable phi by a least squares method,     !
!   without refershing the buffers.                                            !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Vof_Type), target :: Vof
  real                    :: phi (-Vof % pnt_grid % n_bnd_cells:  &
                                   Vof % pnt_grid % n_cells)
  integer, intent(in)     :: i
  real,    intent(out)    :: phii(-Vof % pnt_grid % n_bnd_cells:  &
                                   Vof % pnt_grid % n_cells)
  real,    intent(in)     :: phif  ! phi at the front
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: grid
  type(Field_Type), pointer :: flow
  integer                   :: s, c1, c2
  real                      :: dphi1, dphi2
  real                      :: dx_c1, dy_c1, dz_c1, dx_c2, dy_c2, dz_c2
!-----------------------------[Local parameters]-------------------------------!
  integer, dimension(3,3), parameter :: MAP = reshape((/ 1, 4, 5,  &
                                                         4, 2, 6,  &
                                                         5, 6, 3 /), shape(MAP))
!==============================================================================!

  ! Take alias
  grid => Vof % pnt_grid
  flow => Vof % pnt_flow

  ! Initialize gradients
  phii(1:grid % n_cells) = 0.

  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    ! Take differences as if front doesn't exist
    dphi1 = phi(c2)-phi(c1)
    dphi2 = phi(c2)-phi(c1)
    dx_c1 = grid % dx(s)
    dy_c1 = grid % dy(s)
    dz_c1 = grid % dz(s)
    dx_c2 = grid % dx(s)
    dy_c2 = grid % dy(s)
    dz_c2 = grid % dz(s)

    ! If face is at the front, reduce the extents of the stencil
    if(any(Vof % Front % face_at_elem(1:2,s) .ne. 0)) then
      dphi1 = phif - phi(c1)
      dphi2 = phi(c2) - phif
      dx_c1 = grid % xs(s) - grid % xc(c1)
      dy_c1 = grid % ys(s) - grid % yc(c1)
      dz_c1 = grid % zs(s) - grid % zc(c1)
      dx_c2 = grid % xc(c2) - grid % xs(s)
      dy_c2 = grid % yc(c2) - grid % ys(s)
      dz_c2 = grid % zc(c2) - grid % zs(s)
    end if

    ! On the boundaries
    if(c2 < 0) then
      if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. SYMMETRY) then
        dphi1 = 0.
      end if

      phii(c1) = phii(c1)                                          &
               + dphi1 * (  flow % grad_c2c(MAP(i,1),c1) * dx_c1   &
                          + flow % grad_c2c(MAP(i,2),c1) * dy_c1   &
                          + flow % grad_c2c(MAP(i,3),c1) * dz_c1)
    end if

    ! Inside the domain
    if(c2 > 0) then
      phii(c1) = phii(c1)                                          &
               + dphi1 * (  flow % grad_c2c(MAP(i,1),c1) * dx_c1   &
                          + flow % grad_c2c(MAP(i,2),c1) * dy_c1   &
                          + flow % grad_c2c(MAP(i,3),c1) * dz_c1)

      phii(c2) = phii(c2)                                          &
               + dphi2 * (  flow % grad_c2c(MAP(i,1),c2) * dx_c2   &
                          + flow % grad_c2c(MAP(i,2),c2) * dy_c2   &
                          + flow % grad_c2c(MAP(i,3),c2) * dz_c2)
    end if
  end do

  end subroutine
