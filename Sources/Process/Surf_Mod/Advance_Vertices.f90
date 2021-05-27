!==============================================================================!
  subroutine Advance_Vertices(Surf, Flow)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Surf_Type),  target :: Surf
  class(Field_Type), target :: Flow
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
  type(Vert_Type), pointer :: Vert(:)
  type(Var_Type),  pointer :: u, v, w
  integer,         pointer :: nv
  integer                  :: ver, c, n_verts_in_buffers
  real                     :: max_dis, rx, ry, rz, r, u_v, v_v, w_v
!==============================================================================!

  ! Take aliases
  nv   => Surf % n_verts
  Vert => Surf % Vert
  Grid => Surf % pnt_grid
  u    => Flow % u
  v    => Flow % v
  w    => Flow % w

  max_dis = 0  ! was used for checking

  !----------------------------------------------!
  !   The strategy outlined here, clusters the   !
  !   vertices at the trailing end of surface    !
  !----------------------------------------------!
  do ver = 1, nv

    c = Vert(ver) % cell

    ! Vector from cell center to vertex
    rx = Vert(ver) % x_n - Grid % xc(c)
    ry = Vert(ver) % y_n - Grid % yc(c)
    rz = Vert(ver) % z_n - Grid % zc(c)

    r = sqrt(rx**2 + ry**2 + rz**2)
    max_dis = max(max_dis, r)

    ! Velocity at the vertex
    u_v = u % n(c) + u % x(c) * rx  &
                   + u % y(c) * ry  &
                   + u % z(c) * rz

    v_v = v % n(c) + v % x(c) * rx  &
                   + v % y(c) * ry  &
                   + v % z(c) * rz

    w_v = w % n(c) + w % x(c) * rx  &
                   + w % y(c) * ry  &
                   + w % z(c) * rz


    ! New vertex position
    Vert(ver) % x_n = Vert(ver) % x_n + u_v * Flow % dt
    Vert(ver) % y_n = Vert(ver) % y_n + v_v * Flow % dt
    Vert(ver) % z_n = Vert(ver) % z_n + w_v * Flow % dt

    ! It might have moved to a new cell
    call Vert(ver) % Find_Nearest_Cell(n_verts_in_buffers)
    call Vert(ver) % Find_Nearest_Node()
  end do

  end subroutine
