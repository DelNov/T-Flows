!==============================================================================!
  subroutine Find_Front_Elements(Front, verbose)
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

  !-------------------------------!
  !   Find elements' neighbours   !
  !-------------------------------!
  do s = 1, ns
    c  = side(s) % c
    d  = side(s) % d
    ea = side(s) % ea
    eb = side(s) % eb

    ! Element a
    if(ea > 0) then

      do i_ver = 1, elem(ea) % nv

        ! Get first and second vertex
        v1 = elem(ea) % v(i_ver)
        if(i_ver < elem(ea) % nv) then
          v2 = elem(ea) % v(i_ver+1)
        else
          v2 = elem(ea) % v(1)
        end if

        ! Check if nodes match
        if(v1 .eq. c .and. v2 .eq. d  .or. &
           v2 .eq. c .and. v1 .eq. d) then
          elem(ea) % ns  = elem(ea) % ns  + 1
          elem(ea) % s(elem(ea) % ns)  = s
          if(eb .gt. 0) then
            elem(ea) % nne = elem(ea) % nne + 1
            elem(ea) % e(elem(ea) % nne) = eb
          end if
        end if
      end do

    end if  ! ea > 0

    ! Element b
    if(eb > 0) then

      do i_ver = 1, elem(eb) % nv

        ! Get first and second vertex
        v1 = elem(eb) % v(i_ver)
        if(i_ver < elem(eb) % nv) then
          v2 = elem(eb) % v(i_ver+1)
        else
          v2 = elem(eb) % v(1)
        end if

        ! Check if nodes match
        if(v1 .eq. c .and. v2 .eq. d  .or. &
           v2 .eq. c .and. v1 .eq. d) then
          elem(eb) % ns  = elem(eb) % ns  + 1
          elem(eb) % s(elem(eb) % ns)  = s
          if(ea .gt. 0) then
            elem(eb) % nne = elem(eb) % nne + 1
            elem(eb) % e(elem(eb) % nne) = ea
          end if
        end if
      end do

    end if  ! eb > 0

  end do

  end subroutine
