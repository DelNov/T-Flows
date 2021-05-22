!==============================================================================!
  subroutine Check_Elements(Front, verbose)
!------------------------------------------------------------------------------!
!   Finds connectivity for sides and elements                                  !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Front_Type), target :: Front
  logical                   :: verbose
!-----------------------------------[Locals]-----------------------------------!
  type(Vert_Type), pointer :: Vert(:)
  type(Side_Type), pointer :: side(:)
  type(Elem_Type), pointer :: elem(:)
  integer,         pointer :: nv, ns, ne
  integer                  :: cnt_one, cnt_two
  integer                  :: e, eb, ea, c, d, c1, c2, d1, d2, s, n_side
  integer                  :: ss, sum_ijk, sum_cd, i_ver, j_ver, k_ver
  integer                  :: i_s, v1, v2, v3
  integer                  :: ne_tot, ns_tot, cnt_one_tot, cnt_two_tot
  integer, allocatable     :: ci(:), di(:), ei(:), ni(:)
!==============================================================================!

  ! Take aliases
  nv   => Front % n_verts
  ns   => Front % n_sides
  ne   => Front % n_elems
  Vert => Front % Vert
  side => Front % side
  elem => Front % elem

  ! Checking
  do e = 1, ne
    sum_ijk = 0
    sum_cd  = 0

    do i_ver = 1, elem(e) % nv
      sum_ijk = sum_ijk + elem(e) % v(i_ver)
    end do

    do i_s = 1, elem(e) % ns
      sum_cd = sum_cd + side(elem(e) % s(i_s)) % c  &
                      + side(elem(e) % s(i_s)) % d
    end do

    if( sum_cd / sum_ijk .ne. 2 ) then
      print *, '# ERROR in forming elements'' neighbours!'
      stop
    end if
  end do

  end subroutine
