!==============================================================================!
  subroutine Front_Mod_Handle_4_Points(front, surf_v)
!------------------------------------------------------------------------------!
!   Surface intersects cell at four points                                     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Front_Type), target :: front
  real                     :: surf_v(3)
!-----------------------------------[Locals]-----------------------------------!
  type(Vert_Type), pointer :: vert(:)
  type(Elem_Type), pointer :: elem(:)
  integer,         pointer :: nv, ne
  integer                  :: ver(4), loop
  real                     :: v_21(3), v_31(3), v_41(3)
  real                     :: tri_p_123(3), tri_p_134(3)
  integer                  :: permutations(4, 6)
!=============================================================================!

  ! Take aliases
  nv   => front % n_verts
  ne   => front % n_elems
  vert => front % vert
  elem => front % elem

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

    tri_p_123 = Math_Mod_Cross_Product(v_21, v_31)
    tri_p_134 = Math_Mod_Cross_Product(v_31, v_41)

    if(dot_product(surf_v, tri_p_123) > 0.0 .and.  &
       dot_product(surf_v, tri_p_134) > 0.0) then
      exit
    end if

  end do

  !--------------------------------------------------------------------!
  !   You've found a permutation which works, store new elements now   !
  !--------------------------------------------------------------------!

  ! One new element with three vertices
  ne = ne + 1
  elem(ne) % v(1) = ver(1)
  elem(ne) % v(2) = ver(2)
  elem(ne) % v(3) = ver(3)

  ! Second new element with three vertices
  ne = ne + 1
  elem(ne) % v(1) = ver(1)
  elem(ne) % v(2) = ver(3)
  elem(ne) % v(3) = ver(4)

  end subroutine
