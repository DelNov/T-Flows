!==============================================================================!
  subroutine Grad_Component_Cells_To_Nodes(Flow, phic, phin, i, phii)
!------------------------------------------------------------------------------!
!   Calculates gradient of generic variable phi by the least squares method.   !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Field_Type), target :: Flow
  real                      :: phic( -Flow % pnt_grid % n_bnd_cells  &
                                     :Flow % pnt_grid % n_cells)
  real                      :: phin(1:Flow % pnt_grid % n_nodes)
  integer, intent(in)       :: i
  real,    intent(out)      :: phii(1:Flow % pnt_grid % n_nodes)
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
  integer                  :: n, c, i_cel
  real                     :: dx, dy, dz, dphi
!-----------------------------[Local parameters]-------------------------------!
  integer, dimension(3,3), parameter :: MAP = reshape((/ 1, 4, 5,  &
                                                         4, 2, 6,  &
                                                         5, 6, 3 /), shape(MAP))
!==============================================================================!

  ! Take alias
  Grid => Flow % pnt_grid

  ! Do not refresh buffers and nodal values because
  ! curvature at the walls are imposed outside

  ! Initialize gradients at nodes
  phii(1:Grid % n_nodes) = 0.

  do n = 1, Grid % n_nodes

    ! Browse through cell's nodes
    do i_cel = 1, Grid % nodes_n_cells(n)

      c  = Grid % nodes_c(i_cel, n)
      dx = Grid % xc(c) - Grid % xn(n)
      dy = Grid % yc(c) - Grid % yn(n)
      dz = Grid % zc(c) - Grid % zn(n)

      dphi = phic(c) - phin(n)

      phii(n) = phii(n)                                      &
              + dphi * (  Flow % grad_c2n(MAP(i,1),n) * dx   &
                        + Flow % grad_c2n(MAP(i,2),n) * dy   &
                        + Flow % grad_c2n(MAP(i,3),n) * dz)
    end do

  end do

  end subroutine
