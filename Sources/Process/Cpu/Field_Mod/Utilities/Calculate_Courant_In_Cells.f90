!==============================================================================!
  subroutine Calculate_Courant_In_Cells(Flow, courant)
!------------------------------------------------------------------------------!
!   Calculates Courant number by comparing the path imaginar particle makes    !
!   in a time step, with the maximum span of the cell dimension in the dir-    !
!   estion of the velocity.  To estimate that, imagine that cell center and    !
!   velocity stemming from it define a plane.  Distance of cell nodes from     !
!   this plane will be used to estimete the dimension of the cell in direc-    !
!   tion of the velocity.                                                      !
!------------------------------------------------------------------------------!
  implicit none
!--------------------------------[Arguments]-----------------------------------!
  class(Field_Type), target :: Flow
  real                      :: courant(- Flow % pnt_grid % n_bnd_cells  &
                                       : Flow % pnt_grid % n_cells)
!----------------------------------[Locals]------------------------------------!
  type(Grid_Type), pointer :: Grid
  type(Var_Type),  pointer :: u, v, w
  integer                  :: c, i_nod, n
  real                     :: u_m, n_x, n_y, n_z  ! normalized velocity
  real                     :: d_x, d_y, d_z       ! normalized velocity
  real                     :: dist_pos, dist_neg
!==============================================================================!

  ! Take aliases
  Grid => Flow % pnt_grid
  call Flow % Alias_Momentum(u, v, w)

  !--------------------------------!
  !   Browse through cells first   !
  !--------------------------------!
  do c = 1, Grid % n_cells

    ! Normalize the velocity
    u_m = norm2( (/u % n(c), v % n(c), w % n(c)/) )
    n_x = u % n(c) / (u_m + TINY)
    n_y = v % n(c) / (u_m + TINY)
    n_z = w % n(c) / (u_m + TINY)

    ! Initialize distance from each side of the velocity
    dist_pos = -HUGE
    dist_neg = +HUGE

    ! Browse through nodes of the cell
    do i_nod = 1, abs(Grid % cells_n_nodes(c))
      n = Grid % cells_n(i_nod, c)

      ! Connection between the node and the cell center
      d_x = Grid % xn(n) - Grid % xc(c)
      d_y = Grid % yn(n) - Grid % yc(c)
      d_z = Grid % zn(n) - Grid % zc(c)

      dist_pos = max(dist_pos, d_x*n_x + d_y*n_y + d_z*n_z)
      dist_neg = min(dist_neg, d_x*n_x + d_y*n_y + d_z*n_z)
    end do

    ! With maximum span of the cell in the direction of velocity
    courant(c) = u_m * Flow % dt / (abs(dist_pos) + abs(dist_neg) + TINY)

  end do
  call Grid % Exchange_Cells_Real(courant(-Grid % n_bnd_cells:Grid % n_cells))

  end subroutine
