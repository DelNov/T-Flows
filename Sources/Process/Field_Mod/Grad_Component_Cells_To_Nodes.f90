!==============================================================================!
  subroutine Field_Mod_Grad_Component_Cells_To_Nodes(flow, phic, phin, i, phii)
!------------------------------------------------------------------------------!
!   Calculates gradient of generic variable phi by the least squares method.   !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type), target :: flow
  real                     :: phic( -flow % pnt_grid % n_bnd_cells  &
                                    :flow % pnt_grid % n_cells)
  real                     :: phin(1:flow % pnt_grid % n_nodes)
  integer, intent(in)      :: i
  real,    intent(out)     :: phii(1:flow % pnt_grid % n_nodes)
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer           :: grid
  integer                            :: n, c, i_cel
  real                               :: dx, dy, dz, dphi
  logical                            :: imp_sym
!-----------------------------[Local parameters]-------------------------------!
  integer, dimension(3,3), parameter :: MAP = reshape((/ 1, 4, 5,  &
                                                         4, 2, 6,  &
                                                         5, 6, 3 /), shape(MAP))
!==============================================================================!

  ! Take alias
  grid => flow % pnt_grid

  ! Do not refresh buffers and nodal values because
  ! curvature at the walls are imposed outside

  ! Initialize gradients at nodes
  phii(1:grid % n_nodes) = 0.

  do n = 1, grid % n_nodes

    ! Browse through cell's nodes
    do i_cel = 1, grid % nodes_n_cells(n)

      c  = grid % nodes_c(i_cel, n)
      dx = grid % xc(c) - grid % xn(n)
      dy = grid % yc(c) - grid % yn(n)
      dz = grid % zc(c) - grid % zn(n)

      dphi = phic(c) - phin(n)

      phii(n) = phii(n)                                      &
              + dphi * (  flow % grad_c2n(MAP(i,1),n) * dx   &
                        + flow % grad_c2n(MAP(i,2),n) * dy   &
                        + flow % grad_c2n(MAP(i,3),n) * dz)
    end do

  end do

  end subroutine
