!==============================================================================!
  subroutine Calculate_Element_Normals(Front)
!------------------------------------------------------------------------------!
!>  This subroutine calculates the normals of elements in a front.
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * Iterating through all elements.                                          !
!   * Initializing the area and surface vectors (sx, sy, sz) of each element   !
!     to zero.                                                                 !
!   * Calculating cross products between vectors formed by vertices of each    !
!     element to determine the area and orientation.                           !
!   * Summing these cross products to find the total surface vector for each   !
!     element.                                                                 !
!   * Dividing the total surface vector by the total area to compute the       !
!     normal vector.                                                           !
!   * Verifying the orientation of each element.                               !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Front_Type), target :: Front  !! parent class
!-----------------------------------[Locals]-----------------------------------!
  type(Vert_Type), pointer :: Vert(:)
  type(Elem_Type), pointer :: Elem(:)
  integer,         pointer :: nv, ne
  integer                  :: e, j, v1, v2, i_ver
  real                     :: surf_v(3)
  real                     :: a(3), b(3), tri_v(3), area_x2
  character(SL)            :: str
!==============================================================================!

  ! Take aliases
  nv   => Front % n_verts
  ne   => Front % n_elems
  Vert => Front % Vert
  Elem => Front % Elem

  !---------------------------------!
  !   Browse through all elements   !
  !---------------------------------!
  do e = 1, ne

    Elem(e) % area = 0.0
    Elem(e) % sx   = 0.0
    Elem(e) % sy   = 0.0
    Elem(e) % sz   = 0.0

    Assert(Elem(e) % nv > 2)

    do i_ver = 1, Elem(e) % nv

      v1 = Elem(e) % v(i_ver)
      if(i_ver < Elem(e) % nv) then
        v2 = Elem(e) % v(i_ver + 1)
      else
        v2 = Elem(e) % v(1)
      end if

      a(1) = Vert(v1) % x_n - Elem(e) % xe
      a(2) = Vert(v1) % y_n - Elem(e) % ye
      a(3) = Vert(v1) % z_n - Elem(e) % ze

      b(1) = Vert(v2) % x_n - Elem(e) % xe
      b(2) = Vert(v2) % y_n - Elem(e) % ye
      b(3) = Vert(v2) % z_n - Elem(e) % ze

      tri_v = Math % Cross_Product(a, b)

      ! Magnitude of the cross product (twice the area)
      area_x2 = sqrt(tri_v(1)**2 + tri_v(2)**2 + tri_v(3)**2)

      Elem(e) % sx = Elem(e) % sx + tri_v(1) * 0.5
      Elem(e) % sy = Elem(e) % sy + tri_v(2) * 0.5
      Elem(e) % sz = Elem(e) % sz + tri_v(3) * 0.5

      Elem(e) % area = Elem(e) % area + area_x2 * 0.5

    end do

    ! Compute normals
    Elem(e) % nx = Elem(e) % sx / Elem(e) % area
    Elem(e) % ny = Elem(e) % sy / Elem(e) % area
    Elem(e) % nz = Elem(e) % sz / Elem(e) % area

    do j = 1, 3

      ! Surface vector
      surf_v(1) = Vert(Elem(e) % v(j)) % smooth_x
      surf_v(2) = Vert(Elem(e) % v(j)) % smooth_y
      surf_v(3) = Vert(Elem(e) % v(j)) % smooth_z

      if(dot_product(surf_v, tri_v) < 0.0) then
        write(str, '(i0.0)')  e
        call Message % Error(60,                                          &
                  ' Element '//trim(str)//' is not properly oriented! ',  &
                  file=__FILE__, line=__LINE__)
      end if
    end do   ! for i, j, k

  end do

  end subroutine
