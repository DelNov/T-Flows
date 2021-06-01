!==============================================================================!
  subroutine Smooth_Surf(Surf, smooth)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Surf_Type), target :: Surf
  type(Var_Type),   target :: smooth
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
  type(Vert_Type), pointer :: Vert(:)
  type(Elem_Type), pointer :: Elem(:)
  type(Side_Type), pointer :: side(:)
  integer                  :: it_smooth, s, v, e, c
  integer                  :: n_verts_in_buffers
  integer,         pointer :: nv, ne, ns
  real                     :: xc, yc, zc
  real                     :: nx, ny, nz, phi_m, dx, dy, dz, dm
  real                     :: phi_v
!==============================================================================!

  ! Take aliases
  nv   => Surf % n_verts
  ne   => Surf % n_elems
  ns   => Surf % n_sides
  Vert => Surf % Vert
  Elem => Surf % Elem
  side => Surf % side
  Grid => Surf % pnt_grid

  !-----------------------------------------------------------!
  !   Count number of neighbouring elements for each vertex   !
  !-----------------------------------------------------------!
  Vert(1:nv) % nne = 0
  do e = 1, ne
    Vert(Elem(e) % v(1)) % nne = Vert(Elem(e) % v(1)) % nne + 1
    Vert(Elem(e) % v(2)) % nne = Vert(Elem(e) % v(2)) % nne + 1
    Vert(Elem(e) % v(3)) % nne = Vert(Elem(e) % v(3)) % nne + 1
  end do

  call Surf % Find_Boundaries()

  !---------------------------------------------------------------!
  !   Initialize sums of x, y and z coordinates for each vertex   !
  !---------------------------------------------------------------!
  Vert(1:nv) % sumx = 0.0
  Vert(1:nv) % sumy = 0.0
  Vert(1:nv) % sumz = 0.0

  !-------------------------------------!
  !   Go through smoothing iterations   !
  !-------------------------------------!
  do it_smooth = 1, 3

    ! Compute sums of all coordinates for vertices
    do s = 1, ns
      Vert(side(s) % c) % sumx =  &
      Vert(side(s) % c) % sumx + Vert(side(s) % d) % x_n
      Vert(side(s) % c) % sumy =  &
      Vert(side(s) % c) % sumy + Vert(side(s) % d) % y_n
      Vert(side(s) % c) % sumz =  &
      Vert(side(s) % c) % sumz + Vert(side(s) % d) % z_n

      Vert(side(s) % d) % sumx =  &
      Vert(side(s) % d) % sumx + Vert(side(s) % c) % x_n
      Vert(side(s) % d) % sumy =  &
      Vert(side(s) % d) % sumy + Vert(side(s) % c) % y_n
      Vert(side(s) % d) % sumz =  &
      Vert(side(s) % d) % sumz + Vert(side(s) % c) % z_n
    end do

    ! Move the vertices to their new positions
    do v = 1, nv

      ! This is how I check if vertex is on a boundary.  Maybe there
      ! is a more spothisticated way to do it, but this works so far
      if( .not. Vert(v) % boundary ) then
        Vert(v) % x_n = Vert(v) % sumx / Vert(v) % nne
        Vert(v) % y_n = Vert(v) % sumy / Vert(v) % nne
        Vert(v) % z_n = Vert(v) % sumz / Vert(v) % nne
        call Surf % Vert(v) % Find_Nearest_Cell(n_verts_in_buffers)
        call Surf % Vert(v) % Find_Nearest_Node()
        c = Vert(v) % cell
        Surf % Vert(v) % smooth   = smooth % n(c)
        Surf % Vert(v) % smooth_x = smooth % x(c)
        Surf % Vert(v) % smooth_y = smooth % y(c)
        Surf % Vert(v) % smooth_z = smooth % z(c)
      end if

    end do

    ! Re-initalize the sums
    Vert(1:nv) % sumx = 0.0
    Vert(1:nv) % sumy = 0.0
    Vert(1:nv) % sumz = 0.0

    ! Correct vertex position
    do v = 1, nv

      ! This is how I check if vertex is on a boundary.  Maybe there
      ! is a more spothisticated way to do it, but this works so far
      if( .not. Vert(v) % boundary ) then

        c = Vert(v) % cell

        ! Cell coordinates
        xc = Grid % xc(c)
        yc = Grid % yc(c)
        zc = Grid % zc(c)

        ! Surface normal
        phi_m = sqrt(  Vert(v) % smooth_x**2  &
                     + Vert(v) % smooth_y**2  &
                     + Vert(v) % smooth_z**2)
        nx = Vert(v) % smooth_x / phi_m
        ny = Vert(v) % smooth_y / phi_m
        nz = Vert(v) % smooth_z / phi_m

        ! Value at current vertex position
        dx = Vert(v) % x_n - xc
        dy = Vert(v) % y_n - yc
        dz = Vert(v) % z_n - zc
        phi_v = Vert(v) % smooth + dx * Vert(v) % smooth_x  &
                                 + dy * Vert(v) % smooth_y  &
                                 + dz * Vert(v) % smooth_z

        dm = (0.5 - phi_v) / (  Vert(v) % smooth_x*nx  &
                              + Vert(v) % smooth_y*ny  &
                              + Vert(v) % smooth_z*nz)

        dx = dm * nx
        dy = dm * ny
        dz = dm * nz

        ! Move vertex in the surface normal direction
        Vert(v) % x_n = Vert(v) % x_n + dx
        Vert(v) % y_n = Vert(v) % y_n + dy
        Vert(v) % z_n = Vert(v) % z_n + dz

        dx = Vert(v) % x_n - xc
        dy = Vert(v) % y_n - yc
        dz = Vert(v) % z_n - zc
        phi_v = Vert(v) % smooth + dx * Vert(v) % smooth_x  &
                                 + dy * Vert(v) % smooth_y  &
                                 + dz * Vert(v) % smooth_z

        call Surf % Vert(v) % Find_Nearest_Cell(n_verts_in_buffers)
        call Surf % Vert(v) % Find_Nearest_Node()
        c = Vert(v) % cell
        Surf % Vert(v) % smooth   = smooth % n(c)
        Surf % Vert(v) % smooth_x = smooth % x(c)
        Surf % Vert(v) % smooth_y = smooth % y(c)
        Surf % Vert(v) % smooth_z = smooth % z(c)

      end if  ! if vertex is on a boundary

    end do

  end do  ! smoothing iterations

  ! Make surface exactly spherical for the case of rising bubble in:
  ! Tests/Vof/Rising_Bubble
  ! do v = 1, nv
  !   dm = norm2( (/Vert(v) % x_n, Vert(v) % y_n, Vert(v) % z_n/) )
  !   Vert(v) % x_n = Vert(v) % x_n * 0.25 / dm
  !   Vert(v) % y_n = Vert(v) % y_n * 0.25 / dm
  !   Vert(v) % z_n = Vert(v) % z_n * 0.25 / dm
  ! end do

  end subroutine
