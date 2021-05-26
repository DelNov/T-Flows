!==============================================================================!
  subroutine Handle_4_Points(Front, surf_v, enforce_triangles)
!------------------------------------------------------------------------------!
!   Surface intersects cell at four points                                     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Front_Type), target :: Front
  real                      :: surf_v(3)
  logical                   :: enforce_triangles
!-----------------------------------[Locals]-----------------------------------!
  type(Vert_Type), pointer :: Vert(:)
  type(Elem_Type), pointer :: Elem(:)
  integer,         pointer :: nv, ne
  integer                  :: ver(4), loop
  real                     :: v_21(3), v_31(3), v_41(3)
  real                     :: tri_p_123(3), tri_p_134(3)
  integer                  :: permutations(4, 6)
!=============================================================================!

  ! Take aliases
  nv   => Front % n_verts
  ne   => Front % n_elems
  Vert => Front % Vert
  Elem => Front % Elem

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

    v_21(1) = Vert(ver(2)) % x_n - Vert(ver(1)) % x_n
    v_21(2) = Vert(ver(2)) % y_n - Vert(ver(1)) % y_n
    v_21(3) = Vert(ver(2)) % z_n - Vert(ver(1)) % z_n

    v_31(1) = Vert(ver(3)) % x_n - Vert(ver(1)) % x_n
    v_31(2) = Vert(ver(3)) % y_n - Vert(ver(1)) % y_n
    v_31(3) = Vert(ver(3)) % z_n - Vert(ver(1)) % z_n

    v_41(1) = Vert(ver(4)) % x_n - Vert(ver(1)) % x_n
    v_41(2) = Vert(ver(4)) % y_n - Vert(ver(1)) % y_n
    v_41(3) = Vert(ver(4)) % z_n - Vert(ver(1)) % z_n

    tri_p_123 = Math % Cross_Product(v_21, v_31)
    tri_p_134 = Math % Cross_Product(v_31, v_41)

    if(dot_product(surf_v, tri_p_123) > 0.0 .and.  &
       dot_product(surf_v, tri_p_134) > 0.0) then
      exit
    end if

  end do

  !--------------------------------------------------------------------!
  !   You've found a permutation which works, store new elements now   !
  !--------------------------------------------------------------------!

  if(.not. enforce_triangles) then

    ! One new element with four vertices
    ne = ne + 1
    Elem(ne) % nv = 4
    Elem(ne) % v(1:4) = ver(1:4)

  else  ! only triangles allows

    ! One new element with three vertices
    ne = ne + 1
    Elem(ne) % nv = 3
    Elem(ne) % v(1) = ver(1)
    Elem(ne) % v(2) = ver(2)
    Elem(ne) % v(3) = ver(3)

    ! Second new element with three vertices
    ne = ne + 1
    Elem(ne) % nv = 3
    Elem(ne) % v(1) = ver(1)
    Elem(ne) % v(2) = ver(3)
    Elem(ne) % v(3) = ver(4)

  end if

  end subroutine
