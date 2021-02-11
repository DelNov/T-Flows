!==============================================================================!
  subroutine Front_Mod_Calculate_Element_Normals(front, phi)
!------------------------------------------------------------------------------!
!   Calculates element normals                                                 !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Front_Type), target :: front
  type(Var_Type),   target :: phi
!-----------------------------------[Locals]-----------------------------------!
  type(Vert_Type), pointer :: vert(:)
  type(Elem_Type), pointer :: elem(:)
  integer,         pointer :: nv, ne
  integer                  :: c, e, j, v1, v2, i_v
  real                     :: surf_v(3)
  real                     :: a(3), b(3), tri_v(3), area_x2
!==============================================================================!

  ! Take aliases
  nv   => front % n_verts
  ne   => front % n_elems
  vert => front % vert
  elem => front % elem

  !---------------------------------!
  !   Browse through all elements   !
  !---------------------------------!
  do e = 1, ne

    elem(e) % area = 0.0
    elem(e) % sx   = 0.0
    elem(e) % sy   = 0.0
    elem(e) % sz   = 0.0

    do i_v = 1, elem(e) % nv

      v1 = elem(e) % v(i_v)
      if(i_v < elem(e) % nv) then
        v2 = elem(e) % v(i_v + 1)
      else
        v2 = elem(e) % v(1)
      end if

      a(1) = vert(v1) % x_n - elem(e) % xe
      a(2) = vert(v1) % y_n - elem(e) % ye
      a(3) = vert(v1) % z_n - elem(e) % ze

      b(1) = vert(v2) % x_n - elem(e) % xe
      b(2) = vert(v2) % y_n - elem(e) % ye
      b(3) = vert(v2) % z_n - elem(e) % ze

      tri_v = Math_Mod_Cross_Product(a, b)

      ! Magnitude of the cross product (twice the area)
      area_x2 = sqrt(tri_v(1)**2 + tri_v(2)**2 + tri_v(3)**2)

      elem(e) % sx = elem(e) % sx + tri_v(1) * 0.5
      elem(e) % sy = elem(e) % sy + tri_v(2) * 0.5
      elem(e) % sz = elem(e) % sz + tri_v(3) * 0.5

      elem(e) % area = elem(e) % area + area_x2 * 0.5

    end do

    ! Compute normals
    elem(e) % nx = elem(e) % sx / elem(e) % area
    elem(e) % ny = elem(e) % sy / elem(e) % area
    elem(e) % nz = elem(e) % sz / elem(e) % area

    do j = 1, 3

      ! Take the closest cell
      c = vert(elem(e) % v(j)) % cell

      ! Surface vector
      surf_v(1) = phi % x(c)
      surf_v(2) = phi % y(c)
      surf_v(3) = phi % z(c)

      if(dot_product(surf_v, tri_v) < 0.0) then
        print *, '# Error, element ', e, 'is not properly oriented!'
        stop
      end if
    end do   ! for i, j, k

  end do

  end subroutine
