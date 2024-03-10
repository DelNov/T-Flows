!==============================================================================!
  subroutine Smooth_Surf(Surf, smooth)
!------------------------------------------------------------------------------!
!>  Smooth_Surf is focused on optimizing the mesh's geometric layout by
!>  smoothing out irregularities and enhancing uniformity. This process
!>  involves adjusting vertex positions based on neighboring vertices and
!>  surface normals, ensuring that the mesh accurately represents the physical
!>  surface without unnecessary undulations or sharp angles.  This functions
!>  is based on TRIPOS's (https://github.com/Niceno/TRIPOS) but since this
!>  version is 3D, couldn't quite be the same.
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * Alias setup:                                                             !
!     - Establishes aliases for pointers to vertices (Vert), elements (Elem),  !
!       sides (side), and the grid (Grid). This setup simplifies code          !
!       navigation and enhances clarity.                                       !
!   * Initialization:                                                          !
!     - Initializes variables for counting neighboring elements (nne) and      !
!       accumulating the sum of coordinates (sumx, sumy, sumz) for each vertex.!
!   * Boundary identification:                                                 !
!     - Calls Find_Boundaries to update information about which vertices are   !
!       on boundaries, crucial for ensuring these vertices are treated         !
!       appropriately during the smoothing process.                            !
!   * Smoothing iterations:                                                    !
!     - Executes multiple iterations (3 in this case) to progressively smooth  !
!       the mesh.                                                              !
!     - For each vertex, computes the sum of coordinates of neighboring        !
!       vertices.                                                              !
!     - Adjusts the position of each vertex based on these sums and the number !
!       of neighboring elements, moving it towards the centroid of its         !
!       neighbors.
!   * Boundary vertex handling:                                                !
!     - Checks if a vertex is on a boundary and excludes it from the smoothing !
!       process to preserve the mesh's external geometry.                      !
!   * Surface normal adjustment:                                               !
!     - After moving vertices, adjusts their positions in the direction of the !
!       surface normal. This step ensures that the vertices align with the     !
!       underlying physical surface.                                           !
!   * Re-distribution of smooth variable:                                      !
!     - After each smoothing iteration, re-distributes the smooth variable,    !
!       which likely represents a physical quantity (e.g., temperature,        !
!       pressure) that needs to be updated according to new vertex positions.  !
!   * Correction of vertex position:                                           !
!     - Corrects the position of each vertex by moving it along the surface    !
!       normal based on the smooth variable. This ensures that the mesh        !
!       conforms to the desired physical characteristics of the surface.       !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Surf_Type), target :: Surf    !! parent class
  type(Var_Type),   target :: smooth  !! smooth variant of the scalar field
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
  type(Vert_Type), pointer :: Vert(:)
  type(Elem_Type), pointer :: Elem(:)
  type(Side_Type), pointer :: side(:)
  integer,         pointer :: nv, ne, ns
  integer                  :: it_smooth, s, v, e, c
  integer                  :: n_verts_in_buffers
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
    n_verts_in_buffers = 0
    do v = 1, nv

      ! This is how I check if vertex is on a boundary.  Maybe there
      ! is a more spothisticated way to do it, but this works so far
      if( .not. Vert(v) % boundary ) then
        Vert(v) % x_n = Vert(v) % sumx / Vert(v) % nne
        Vert(v) % y_n = Vert(v) % sumy / Vert(v) % nne
        Vert(v) % z_n = Vert(v) % sumz / Vert(v) % nne
        call Surf % Vert(v) % Find_Nearest_Cell(n_verts_in_buffers)
        call Surf % Vert(v) % Find_Nearest_Node()
      end if

    end do

    ! Re-distribute "smooth" after moving the vertices
    call Surf % Distribute_Smooth(smooth)
    call Surf % Distribute_Cell_Coords()

    ! Re-initalize the sums
    Vert(1:nv) % sumx = 0.0
    Vert(1:nv) % sumy = 0.0
    Vert(1:nv) % sumz = 0.0

    ! Correct vertex position
    n_verts_in_buffers = 0
    do v = 1, nv

      ! This is how I check if vertex is on a boundary.  Maybe there
      ! is a more spothisticated way to do it, but this works so far
      if( .not. Vert(v) % boundary ) then

        c = Vert(v) % cell

        ! Cell coordinates
        xc = Vert(v) % cell_x
        yc = Vert(v) % cell_y
        zc = Vert(v) % cell_z

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

      end if  ! if vertex is on a boundary

    end do

    ! Re-distribute "smooth" after moving the vertices
    call Surf % Distribute_Smooth(smooth)
    call Surf % Distribute_Cell_Coords()

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
