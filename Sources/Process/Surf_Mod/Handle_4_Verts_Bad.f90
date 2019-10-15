!==============================================================================!
  subroutine Surf_Mod_Handle_4_Verts(surf, surf_v)
!------------------------------------------------------------------------------!
!   Places surface where variable phi has value val                            !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Surf_Type), target :: surf
  real                    :: surf_v(3)
!-----------------------------------[Locals]-----------------------------------!
  integer :: ver(4)
  real    :: a(3), b(3), c(3), d(3), tri_v_ad(3), tri_v_cb(3)
!==============================================================================!

  ! Take aliases
  nv   => surf % n_verts
  ne   => surf % n_elems
  vert => surf % vert
  elem => surf % elem

  ! Make the first assumption
  ver(1) = surf % n_verts - 3
  ver(2) = surf % n_verts - 2
  ver(3) = surf % n_verts - 1
  ver(4) = surf % n_verts

  !-------------------------------------------!
  !   Try to find a permutation which works   !
  !-------------------------------------------!
  do

    !        4-----c-----3
    !       /           /
    !      d           b
    !     /           /
    !    1-----a-----2

    ! a x d and c x b
    a(1) = surf % vert(ver(2)) % x_n - surf % vert(ver(1)) % x_n
    a(2) = surf % vert(ver(2)) % y_n - surf % vert(ver(1)) % y_n
    a(3) = surf % vert(ver(2)) % z_n - surf % vert(ver(1)) % z_n

    d(1) = surf % vert(ver(4)) % x_n - surf % vert(ver(1)) % x_n
    d(2) = surf % vert(ver(4)) % y_n - surf % vert(ver(1)) % y_n
    d(3) = surf % vert(ver(4)) % z_n - surf % vert(ver(1)) % z_n

    c(1) = surf % vert(ver(4)) % x_n - surf % vert(ver(3)) % x_n
    c(2) = surf % vert(ver(4)) % y_n - surf % vert(ver(3)) % y_n
    c(3) = surf % vert(ver(4)) % z_n - surf % vert(ver(3)) % z_n

    b(1) = surf % vert(ver(2)) % x_n - surf % vert(ver(3)) % x_n
    b(2) = surf % vert(ver(2)) % y_n - surf % vert(ver(3)) % y_n
    b(3) = surf % vert(ver(2)) % z_n - surf % vert(ver(3)) % z_n

    tri_v_ad = Math_Mod_Cross(a, d)
    tri_v_cb = Math_Mod_Cross(c, b)

    ! Order is good (signs are the same, but orientation is bad
    if(dot_product(surf_v, tri_v_ad) < 0.0 .and.  &
       dot_product(surf_v, tri_v_cb) < 0.0) then
      call Swap_Int(ver(2), ver(4))
    end if

    if(dot_product(surf_v, tri_v_ad) > 0.0 .and.  &
       dot_product(surf_v, tri_v_cb) < 0.0) then
      call Swap_Int(ver(3), ver(4))
    end if

    if(dot_product(surf_v, tri_v_ad) < 0.0 .and.  &
       dot_product(surf_v, tri_v_cb) > 0.0) then
      call Swap_Int(ver(3), ver(4))
    end if

    if(dot_product(surf_v, tri_v_ad) > 0.0 .and.  &
       dot_product(surf_v, tri_v_cb) > 0.0) then
      exit
    end if

  end do

  !--------------------------------------------------------------------!
  !   You've found a permutation which works, store new elements now   !
  !--------------------------------------------------------------------!

  ! One new element with three vertices
  surf % n_elems = surf % n_elems + 1
  surf % elem(surf % n_elems) % i = ver(1)
  surf % elem(surf % n_elems) % j = ver(2)
  surf % elem(surf % n_elems) % k = ver(3)

  ! Second new element with three vertices
  surf % n_elems = surf % n_elems + 1
  surf % elem(surf % n_elems) % i = ver(1)
  surf % elem(surf % n_elems) % j = ver(3)
  surf % elem(surf % n_elems) % k = ver(4)

  end subroutine
