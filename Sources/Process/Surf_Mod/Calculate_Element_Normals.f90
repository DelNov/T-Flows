!==============================================================================!
  subroutine Surf_Mod_Calculate_Element_Normals(surf, phi)
!------------------------------------------------------------------------------!
!   Calculates element normals                                                 !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Surf_Type), target :: surf
  type(Var_Type),  target :: phi
!-----------------------------------[Locals]-----------------------------------!
  type(Vert_Type), pointer :: vert(:)
  type(Elem_Type), pointer :: elem(:)
  integer,         pointer :: nv, ne
  integer                  :: c, e, j
  real                     :: surf_v(3)
  real                     :: a(3), b(3), tri_v(3), area_x2
!==============================================================================!

  ! Take aliases
  nv   => surf % n_verts
  ne   => surf % n_elems
  vert => surf % vert
  elem => surf % elem

  !---------------------------------!
  !   Browse through all elements   !   (in the future, only for this_proc)
  !---------------------------------!
  do e = 1, ne

    a(1) = vert(elem(e) % v(2)) % x_n - vert(elem(e) % v(1)) % x_n
    a(2) = vert(elem(e) % v(2)) % y_n - vert(elem(e) % v(1)) % y_n
    a(3) = vert(elem(e) % v(2)) % z_n - vert(elem(e) % v(1)) % z_n

    b(1) = vert(elem(e) % v(3)) % x_n - vert(elem(e) % v(1)) % x_n
    b(2) = vert(elem(e) % v(3)) % y_n - vert(elem(e) % v(1)) % y_n
    b(3) = vert(elem(e) % v(3)) % z_n - vert(elem(e) % v(1)) % z_n

    tri_v = Math_Mod_Cross_Product(a, b)

    ! Magnitude of the cross product (twice the area)
    area_x2 = sqrt(tri_v(1)**2 + tri_v(2)**2 + tri_v(3)**2)

    elem(e) % nx = tri_v(1) / area_x2
    elem(e) % ny = tri_v(2) / area_x2
    elem(e) % nz = tri_v(3) / area_x2

    do j = 1, 3

      ! Take the closest cell
      if(j .eq. 1) c = vert(elem(e) % v(1)) % cell
      if(j .eq. 2) c = vert(elem(e) % v(2)) % cell
      if(j .eq. 3) c = vert(elem(e) % v(3)) % cell

      ! Surface vector
      surf_v(1) = phi % x(c)
      surf_v(2) = phi % y(c)
      surf_v(3) = phi % z(c)

      if(dot_product(surf_v, tri_v) < 0.0) then
        print *, '# Error, element ', e, 'is not properly oriented!'
!       stop
      end if
    end do   ! for i, j, k

  end do

  end subroutine
