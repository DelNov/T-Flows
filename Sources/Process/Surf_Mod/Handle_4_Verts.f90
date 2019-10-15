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
  type(Vert_Type), pointer :: vert(:)
  type(Elem_Type), pointer :: elem(:)
  integer,         pointer :: nv, ne
  integer                  :: ver(4), loop
  real                     :: v_21(3), v_31(3), v_41(3)
  real                     :: tri_v_123(3), tri_v_134(3)
  integer                  :: permutations(4, 6)
!=============================================================================!

  ! Take aliases
  nv   => surf % n_verts
  ne   => surf % n_elems
  vert => surf % vert
  elem => surf % elem

  permutations = reshape((/ &
    0, 1, 2, 3,  &
    0, 2, 1, 3,  &
    0, 3, 1, 2,  &
    0, 1, 3, 2,  &
    0, 2, 3, 1,  &
    0, 3, 2, 1   &
  /), shape(permutations))

  !-------------------------------------------!
  !   Try to find a permutation which works   !
  !-------------------------------------------!
  do loop = 1, 6

    ver(1) = nv - permutations(1, loop)
    ver(2) = nv - permutations(2, loop)
    ver(3) = nv - permutations(3, loop)
    ver(4) = nv - permutations(4, loop)

    v_21(1) = vert(ver(2)) % x_n - vert(ver(1)) % x_n
    v_21(2) = vert(ver(2)) % y_n - vert(ver(1)) % y_n
    v_21(3) = vert(ver(2)) % z_n - vert(ver(1)) % z_n

    v_31(1) = vert(ver(3)) % x_n - vert(ver(1)) % x_n
    v_31(2) = vert(ver(3)) % y_n - vert(ver(1)) % y_n
    v_31(3) = vert(ver(3)) % z_n - vert(ver(1)) % z_n

    v_41(1) = vert(ver(4)) % x_n - vert(ver(1)) % x_n
    v_41(2) = vert(ver(4)) % y_n - vert(ver(1)) % y_n
    v_41(3) = vert(ver(4)) % z_n - vert(ver(1)) % z_n

    tri_v_123 = Math_Mod_Cross_Product(v_21, v_31)
    tri_v_134 = Math_Mod_Cross_Product(v_31, v_41)

    if(dot_product(surf_v, tri_v_123) > 0.0 .and.  &
       dot_product(surf_v, tri_v_134) > 0.0) then
      exit
    end if

  end do

  !--------------------------------------------------------------------!
  !   You've found a permutation which works, store new elements now   !
  !--------------------------------------------------------------------!

  ! One new element with three vertices
  ne = ne + 1
  elem(ne) % i = ver(1)
  elem(ne) % j = ver(2)
  elem(ne) % k = ver(3)

  ! Second new element with three vertices
  ne = ne + 1
  elem(ne) % i = ver(1)
  elem(ne) % j = ver(3)
  elem(ne) % k = ver(4)

  end subroutine
