!==============================================================================!
  subroutine Field_Mod_Grad_Component_Nodes_To_Cells(flow, phic, phin, i, phii)
!------------------------------------------------------------------------------!
!   Calculates gradient of generic variable phi by a least squares method.     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type), target :: flow
  real                     :: phic( -flow % pnt_grid % n_bnd_cells  &
                                    :flow % pnt_grid % n_cells)
  real                     :: phin(1:flow % pnt_grid % n_nodes)
  integer                  :: i
  real,   intent(out)      :: phii( -flow % pnt_grid % n_bnd_cells  &
                                    :flow % pnt_grid % n_cells)
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer           :: grid
  integer                            :: n, c, i_nod
  real                               :: dx, dy, dz, dphi
  logical                            :: imp_sym
!-----------------------------[Local parameters]-------------------------------!
  integer, dimension(3,3), parameter :: MAP = reshape((/ 1, 4, 5,  &
                                                         4, 2, 6,  &
                                                         5, 6, 3 /), shape(MAP))
!==============================================================================!

  ! Take alias
  grid => flow % pnt_grid

  ! Refresh buffers and nodal values
  call Grid_Mod_Exchange_Cells_Real(grid, phic)
  call Field_Mod_Interpolate_Nodes_To_Cells(flow, phin, phic)

  ! Initialize gradients
  phii(1:grid % n_cells) = 0.

  do c = 1, grid % n_cells

    ! Browse through cell's nodes
    do i_nod = 1, grid % cells_n_nodes(c)

     n  = grid % cells_n(i_nod, c)
     dx = grid % xn(n) - grid % xc(c)
     dy = grid % yn(n) - grid % yc(c)
     dz = grid % zn(n) - grid % zc(c)

     dphi = phin(n) - phic(c)

     phii(c) = phii(c)                                      &
             + dphi * (  flow % grad_n2c(MAP(i,1),c) * dx   &
                       + flow % grad_n2c(MAP(i,2),c) * dy   &
                       + flow % grad_n2c(MAP(i,3),c) * dz)
    end do

  end do

  call Grid_Mod_Exchange_Cells_Real(grid, phii)

  end subroutine
