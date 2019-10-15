!==============================================================================!
  subroutine Surf_Mod_Smooth(surf, phi, val)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Surf_Type), target :: surf
  type(Var_Type),  target :: phi
  real                    :: val
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: grid
  type(Vert_Type), pointer :: vert(:)
  type(Elem_Type), pointer :: elem(:)
  type(Side_Type), pointer :: side(:)
  integer                  :: it_smooth, it_corr, s, v, e, c
  integer                  :: n_verts_in_buffers
  integer,         pointer :: nv, ne, ns
  real                     :: xc, yc, zc
  real                     :: phi_x, phi_y, phi_z, phi_m, dx, dy, dz, dm
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

  print *, 'initial vertex position:         ', vert(1) % x_n, vert(1) % y_n, vert(1) % z_n

  !---------------------------------------------------------!
  !   Count number of neighbouring elements for each node   !
  !---------------------------------------------------------!
  vert(1:nv) % nne = 0
  do e = 1, ne
    vert(elem(e) % i) % nne = vert(elem(e) % i) % nne + 1
    vert(elem(e) % j) % nne = vert(elem(e) % j) % nne + 1
    vert(elem(e) % k) % nne = vert(elem(e) % k) % nne + 1
  end do

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
      vert(v) % x_n = vert(v) % sumx / vert(v) % nne
      vert(v) % y_n = vert(v) % sumy / vert(v) % nne
      vert(v) % z_n = vert(v) % sumz / vert(v) % nne
      call Surf_Mod_Find_Nearest_Cell(surf, v, n_verts_in_buffers)
      call Surf_Mod_Find_Nearest_Node(surf, v)
    end do
    vert(1:nv) % sumx = 0.0
    vert(1:nv) % sumy = 0.0
    vert(1:nv) % sumz = 0.0

    ! Correct vertex position
    do v = 1, nv 
      c = vert(v) % cell

      ! Cell coordinates
      xc = grid % xc(c)
      yc = grid % yc(c)
      zc = grid % zc(c)

      ! Gradient's unit vector
      phi_m = 1.0 ! sqrt(phi % x(c)**2 + phi % y(c)**2 + phi % z(c)**2)
      phi_x = phi % x(c) / phi_m
      phi_y = phi % y(c) / phi_m
      phi_z = phi % z(c) / phi_m

      do it_corr = 1, 3

        ! Vertex coordinates
        dx = vert(v) % x_n - xc
        dy = vert(v) % y_n - yc
        dz = vert(v) % z_n - zc

        val_v = phi % n(c) + dx * phi_x  &
                           + dy * phi_y  &
                           + dz * phi_z
        print *, val_v

        ! Find out how much it_smooth should move to have correct value
        dx = phi_x * (val - val_v)
        dy = phi_y * (val - val_v)
        dz = phi_z * (val - val_v)

        vert(v) % x_n = vert(v) % x_n + dx
        vert(v) % y_n = vert(v) % y_n + dy
        vert(v) % z_n = vert(v) % z_n + dz

!#      call Surf_Mod_Find_Nearest_Cell(surf, v, n_verts_in_buffers)
!#      call Surf_Mod_Find_Nearest_Node(surf, v)
      end do

    end do

  end do  ! smoothing iterations

  end subroutine
