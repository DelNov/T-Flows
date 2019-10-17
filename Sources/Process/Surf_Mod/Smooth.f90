!==============================================================================!
  subroutine Surf_Mod_Smooth(surf, phi, val_e)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Surf_Type), target :: surf
  type(Var_Type),  target :: phi
  real                    :: val_e
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: grid
  type(Vert_Type), pointer :: vert(:)
  type(Elem_Type), pointer :: elem(:)
  type(Side_Type), pointer :: side(:)
  integer                  :: it_smooth, s, v, e, c
  integer                  :: n_verts_in_buffers
  integer,         pointer :: nv, ne, ns
  real                     :: xc, yc, zc
  real                     :: nx, ny, nz, phi_m, dx, dy, dz, dm
  real                     :: val_v
!==============================================================================!

  ! Take aliases
  nv   => surf % n_verts
  ne   => surf % n_elems
  ns   => surf % n_sides
  vert => surf % vert
  elem => surf % elem
  side => surf % side
  grid => surf % pnt_grid

  !---------------------------------------------------------!
  !   Count number of neighbouring elements for each node   !
  !---------------------------------------------------------!
  vert(1:nv) % nne = 0
  do e = 1, ne
    vert(elem(e) % i) % nne = vert(elem(e) % i) % nne + 1
    vert(elem(e) % j) % nne = vert(elem(e) % j) % nne + 1
    vert(elem(e) % k) % nne = vert(elem(e) % k) % nne + 1
  end do

  call Surf_Mod_Find_Boundaries(surf)

  !---------------------------------------------------------------!
  !   Initialize sums of x, y and z coordinates for each vertex   !
  !---------------------------------------------------------------!
  vert(1:nv) % sumx = 0.0
  vert(1:nv) % sumy = 0.0
  vert(1:nv) % sumz = 0.0

  !-------------------------------------!
  !   Go through smoothing iterations   !
  !-------------------------------------!
  do it_smooth = 1, 3

    ! Compute sums of all coordinates for nodes
    do s = 1, ns
      vert(side(s) % c) % sumx =  &
      vert(side(s) % c) % sumx + vert(side(s) % d) % x_n
      vert(side(s) % c) % sumy =  &
      vert(side(s) % c) % sumy + vert(side(s) % d) % y_n
      vert(side(s) % c) % sumz =  &
      vert(side(s) % c) % sumz + vert(side(s) % d) % z_n

      vert(side(s) % d) % sumx =  &
      vert(side(s) % d) % sumx + vert(side(s) % c) % x_n
      vert(side(s) % d) % sumy =  &
      vert(side(s) % d) % sumy + vert(side(s) % c) % y_n
      vert(side(s) % d) % sumz =  &
      vert(side(s) % d) % sumz + vert(side(s) % c) % z_n
    end do

    ! Move the vertices to their new positions
    do v = 1, nv

      ! This is how I check if vertex is on a boundary.  Maybe there
      ! is a more spothisticated way to do it, but this works so far
      if( .not. vert(v) % boundary ) then
        vert(v) % x_n = vert(v) % sumx / vert(v) % nne
        vert(v) % y_n = vert(v) % sumy / vert(v) % nne
        vert(v) % z_n = vert(v) % sumz / vert(v) % nne
        call Surf_Mod_Find_Nearest_Cell(surf, v, n_verts_in_buffers)
        call Surf_Mod_Find_Nearest_Node(surf, v)
      end if

    end do

    ! Re-initalize the sums
    vert(1:nv) % sumx = 0.0
    vert(1:nv) % sumy = 0.0
    vert(1:nv) % sumz = 0.0

    ! Correct vertex position
    do v = 1, nv

      ! This is how I check if vertex is on a boundary.  Maybe there
      ! is a more spothisticated way to do it, but this works so far
      if( .not. vert(v) % boundary ) then

        c = vert(v) % cell

        ! Cell coordinates
        xc = grid % xc(c)
        yc = grid % yc(c)
        zc = grid % zc(c)

        ! Surface normal
        phi_m = sqrt(phi % x(c)**2 + phi % y(c)**2 + phi % z(c)**2)
        nx = phi % x(c) / phi_m
        ny = phi % y(c) / phi_m
        nz = phi % z(c) / phi_m

        ! Value at current vertex position
        dx = vert(v) % x_n - xc
        dy = vert(v) % y_n - yc
        dz = vert(v) % z_n - zc
        val_v = phi % n(c) + dx * phi % x(c)  &
                           + dy * phi % y(c)  &
                           + dz * phi % z(c)

        dm = (val_e - val_v) / (phi % x(c)*nx + phi % y(c)*ny + phi % z(c)*nz)

        dx = dm * nx
        dy = dm * ny
        dz = dm * nz

        ! Move vertex in the surface normal direction
        vert(v) % x_n = vert(v) % x_n + dx
        vert(v) % y_n = vert(v) % y_n + dy
        vert(v) % z_n = vert(v) % z_n + dz

        dx = vert(v) % x_n - xc
        dy = vert(v) % y_n - yc
        dz = vert(v) % z_n - zc
        val_v = phi % n(c) + dx * phi % x(c)  &
                           + dy * phi % y(c)  &
                           + dz * phi % z(c)

        call Surf_Mod_Find_Nearest_Cell(surf, v, n_verts_in_buffers)
        call Surf_Mod_Find_Nearest_Node(surf, v)

      end if  ! if vertex is on a boundary

    end do

  end do  ! smoothing iterations

  end subroutine
