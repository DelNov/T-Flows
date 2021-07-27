!==============================================================================!
  subroutine Grad_Component_Nodes_To_Cells(Flow, phic, phin, i, phii)
!------------------------------------------------------------------------------!
!   Calculates gradient of generic variable phi by a least squares method.     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Field_Type), target :: Flow
  real                      :: phic( -Flow % pnt_grid % n_bnd_cells  &
                                     :Flow % pnt_grid % n_cells)
  real                      :: phin(1:Flow % pnt_grid % n_nodes)
  integer, intent(in)       :: i
  real,    intent(out)      :: phii( -Flow % pnt_grid % n_bnd_cells  &
                                     :Flow % pnt_grid % n_cells)
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
  integer                  :: n, c, i_nod
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

  ! Initialize gradients
  phii(1:Grid % n_cells) = 0.

  do c = 1, Grid % n_cells

    ! Browse through cell's nodes
    do i_nod = 1, abs(Grid % cells_n_nodes(c))

     n  = Grid % cells_n(i_nod, c)
     dx = Grid % xn(n) - Grid % xc(c)
     dy = Grid % yn(n) - Grid % yc(c)
     dz = Grid % zn(n) - Grid % zc(c)

     dphi = phin(n) - phic(c)

     phii(c) = phii(c)                                      &
             + dphi * (  Flow % grad_n2c(MAP(i,1),c) * dx   &
                       + Flow % grad_n2c(MAP(i,2),c) * dy   &
                       + Flow % grad_n2c(MAP(i,3),c) * dz)
    end do

  end do

  call Grid_Mod_Exchange_Cells_Real(Grid, phii)

  end subroutine
